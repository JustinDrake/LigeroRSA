#pragma once
#include <vector>
#include <gmpxx.h>

#include "protocol/WeightedUnionFind.hpp"
#include "protocol/ConstraintSystem.hpp"
#include "SecretSharingNTT.hpp"

#include "FiniteFields.hpp"
#include <nfl/params.hpp>

// #define NDEBUG

// this is consistent with max bound size (rNoise)
#define MAX_BITS_SIZE 750

struct intraBlocReferences {
    std::vector<size_t> blocs;
    std::vector<size_t> cols;
};

// Enumerating message types
enum class Statements : uint64_t {
    Generic_Addition,
    Generic_Multiplication,
    Generic_IntraBloc,
    Generic_InterblocBoundingVariables,
    Generic_IntrablocBoundingVariables,
    Generic_PolyEval,
    RSARound3_keyGen,
    RSARound4_preSieving,
    RSARound5_preSieving,
    RSARound6_partialDecrypt,
    RSARound7_candidateGeneration,
    RSARound8_beaversTriples,
    RSARound10_jacobiTest,
    RSARound11_12_jacobiGCD,
    RSACeremony,
    RSACeremony_Rounds3to6,
    RSACeremony_Connecting_Proofs,
    RSACeremony_Bounding_Variables,
    RSACeremony_Equating_Variables,
    RSACeremony_EarlyCommitments,
    GenerateData
};

constexpr int NbPrimesQ         = 21;
constexpr int alphaCanLength    = 6;
constexpr int alphaPSLength     = 6;
constexpr int alphaGCDLength    = 18;

// This factor drives the amount by which we oversample evaluation points
// Indeed, we avoid picking points that are roots of unity; hence there is a non-null
// probability that the evaluation point we pick is inadequate and therefore 
// we need excess randomness
constexpr uint64_t evaluationPointsOverSamplingFactor = 10;

std::vector<unsigned char> hashPublicData() {
    std::vector<unsigned char> tmp(iSeedSize);
    return tmp;
}

size_t determineCSrandomnessSize(Statements expressNPStatement, size_t numberOfEvalPoints) {
    size_t size = 0;
    switch(expressNPStatement) {
        case Statements::RSARound3_keyGen:
        case Statements::RSARound4_preSieving:
        case Statements::RSARound5_preSieving:
        case Statements::RSARound6_partialDecrypt:
        size = numberOfEvalPoints;
        break;

        case Statements::RSACeremony:
        case Statements::RSACeremony_Rounds3to6:
        size = 4 * numberOfEvalPoints;
        break;

        case Statements::RSACeremony_Connecting_Proofs:
        size = 2 * numberOfEvalPoints;
        break;

        case Statements::RSACeremony_EarlyCommitments:
        case Statements::RSACeremony_Equating_Variables:
        case Statements::RSACeremony_Bounding_Variables:
        case Statements::RSARound7_candidateGeneration:
        case Statements::RSARound8_beaversTriples:
        case Statements::RSARound10_jacobiTest:
        case Statements::RSARound11_12_jacobiGCD:
        default:
        break;
    }

    return (size*evaluationPointsOverSamplingFactor);
}

template<typename FieldT>
class expressNPStatement {
public:
    std::vector<size_t>     _siCoefsBlocIds;
    std::vector<size_t>     _eiCoefsBlocIds;
    std::vector<size_t>     _xiCoefsBlocIds;

protected:
    size_t                  _blockSize;
    builder<FieldT>&        _builder;
    PublicData&             _publicData;
    SecretData&             _secretData;
    vars<FieldT>            _variables;
    uint64_t*               _phis;
    std::vector<size_t>     _xprimeIds;
    std::vector<size_t>     _yprimeIds;
    std::vector<size_t>     _zprimeIds;
    intraBlocReferences     _xsharesPSidxs;
    intraBlocReferences     _xsharesPSidxs_2;
    intraBlocReferences     _xsharesPSidxs_3;
    intraBlocReferences     _xsharesPSidxs_4;
    intraBlocReferences     _xsharesCANidxs;
    intraBlocReferences     _xsharesGCDidxs;
    intraBlocReferences     _ysharesPSidxs;
    intraBlocReferences     _ysharesPSidxs_2;
    intraBlocReferences     _ysharesPSidxs_3;
    intraBlocReferences     _ysharesPSidxs_4;
    intraBlocReferences     _ysharesCANidxs;
    intraBlocReferences     _ysharesGCDidxs;
    intraBlocReferences     _zsharesCANidxs;
    uint64_t*&               _randomness;
    std::vector<uint64_t>   _qdivpPS, _qdivpCAN, _qdivpGCD;
    bool                    _csPickRandom;

    auto rephraseAs2Nct(uint64_t d) {
        
        // Rephrase the problem using the nearest 2^n 
        auto [ D, log2D ] = ligero::roundToUpperPowerTwo(2*d);
        FieldT C = FieldT(D-d);

        // defining the constraint
        if (D == 0) {
            if (d == 0) {
                DBG("Limit set to 0!");
            } else {
                DBG("Overflow");
            }

            assert(d > 0);
        } else {
            DBG("log2D = " << log2D);
            DBG("d = " << d);
            DBG("C = " << C.getValue());
            DBG("D = " << D);
            DBG("d + C = D");
        }

        return std::tuple<FieldT, FieldT, FieldT, uint64_t>(d,C,D,log2D);
    }

    std::vector<size_t> addSplitConstant(std::vector<FieldT>& ct) {
        std::vector<size_t> blocksIdx;

        for (size_t idx = 0; idx < ct.size()/this->_blockSize; idx++) {
            blocksIdx.push_back(_builder.add(new constant_constraint<FieldT>(vector<FieldT>(ct.begin() + idx*this->_blockSize, ct.begin() + (idx+1)*this->_blockSize))));
        }

        return blocksIdx;
    }

    size_t addPartialConstant(std::vector<FieldT>& ct, std::vector<size_t>& positions) {

        std::vector<FieldT> partial(this->_blockSize, FieldT(0));
        size_t idx = 0;
        for (auto pos: positions) partial[idx++] = ct[pos];

        return _builder.add(new constant_constraint<FieldT>(partial));
    }

    std::vector<size_t> addSplitVariable(std::vector<FieldT>& var) {
        std::vector<size_t> blocksIdx;

        for (size_t idx = 0; idx < static_cast<size_t>(ceil((float)var.size()/(float)this->_blockSize)); idx++) {
            blocksIdx.push_back(_builder.add(new variable_constraint<FieldT>()));

            if ((idx+1)*this->_blockSize <= var.size()) {
                _variables.add(blocksIdx.back(),vector<FieldT>(var.begin() + idx*this->_blockSize, var.begin() + (idx+1)*this->_blockSize));
            } else {
                std::vector<FieldT> temp(var.begin() + idx*this->_blockSize, var.end());
                temp.resize(this->_blockSize, FieldT(0)); 
                _variables.add(blocksIdx.back(),temp);
            }
        }

        return blocksIdx;
    }

    size_t addIntrablocEquality(size_t blocZeroIdx, std::vector<size_t> blocs_left, std::vector<size_t> cols_left, std::vector<size_t> blocs_right, std::vector<size_t> cols_right) {

        size_t blocCol = 0;
        transformation_constraint<FieldT> *tc = new transformation_constraint<FieldT>;

        for (size_t idx =0; idx < blocs_left.size(); idx++) {

            tc->block.emplace_back(blocs_left[idx]);
            tc->source_position.emplace_back(cols_left[idx]);
            tc->scalar.emplace_back(1);
            tc->target_position.emplace_back(blocCol);

            tc->block.emplace_back(blocs_right[idx]);
            tc->source_position.emplace_back(cols_right[idx]);
            tc->scalar.emplace_back(FieldT(FieldT::getModulus() - 1));
            tc->target_position.emplace_back(blocCol);

            blocCol++;
        }

        // padBloc(tc, blocCol, blocZeroIdx);
        return (_builder.add(tc));
    }

    std::vector<FieldT> embed(std::vector<uint64_t> rawData, size_t offset, size_t size) {
        std::vector<FieldT> tmp;
        tmp.reserve(size);
        for (size_t idx = 0; idx < size; idx ++) {
            tmp.push_back(rawData[idx+offset]);
        }
        return tmp;
    }

    void Preprocessing() {
		    auto ntt = new NTT(this->_publicData.modulusIdx, 65536ull, permut<65536>::compute);

            this->_phis = (uint64_t*)calloc(65536, sizeof(uint64_t));
            this->_phis[1] = 1;
            ntt->ntt(this->_phis);

            mpz_class q_over_p = 1;
            for (auto i = 9; i < 21; i++) {
                q_over_p *= nfl::params<uint64_t>::P[i];
            }

            nfl::poly_p<uint64_t, 65536, 21> qdivp{q_over_p};
            qdivp.ntt_pow_phi();

            //PS
            for (size_t idx = 0; idx < 6; idx++) {
                this->_qdivpPS.emplace_back(qdivp.poly_obj().data()[_publicData.modulusIdx*_publicData.degree + _publicData.indicesPS[idx]]);
            }

            //CAN
            for (size_t idx = 0; idx < 6; idx++) {
                this->_qdivpCAN.emplace_back(qdivp.poly_obj().data()[_publicData.modulusIdx*_publicData.degree + _publicData.indicesCAN[idx]]);
            }

            //GCD
            for (size_t idx = 0; idx < 18; idx++) {
                this->_qdivpGCD.emplace_back(qdivp.poly_obj().data()[_publicData.modulusIdx*_publicData.degree + _publicData.indicesGCD[idx]]);
            }

            delete ntt;
    }

    void NP_RSACeremony_Bounding_Variables(size_t blocZeroIdx, size_t blocOneIdx) {
            std::vector<FieldT> varsXs;
            std::vector<FieldT> varsXplusCs;
            std::vector<FieldT> varsXplusDs;    
            std::vector<FieldT> Cs;
            std::vector<FieldT> ds;

            std::vector<std::vector<FieldT>> bitDecomposition;
            std::vector<size_t> log2Ds;

            size_t sizeIndices = (_publicData.indicesPS.size()+_publicData.indicesCAN.size()+_publicData.indicesGCD.size());

            log2Ds.assign(_publicData.log2Ds.begin(), _publicData.log2Ds.end());
            
            for (size_t idx = 0; idx < sizeIndices; idx++) {
                int ref = _publicData.modulusIdx + idx*NbPrimesQ;
                Cs.emplace_back(_publicData.Cs[ref]);
                ds.emplace_back(_publicData.Ds[ref]);
            }

            if (_secretData.isAvailable) {
                for (size_t idx = 0; idx < sizeIndices; idx++) {
                    int ref = _publicData.modulusIdx + idx*NbPrimesQ;
                    
                    varsXs.emplace_back(_secretData.Xs[ref]);
                    varsXplusCs.emplace_back(_secretData.XplusCs[ref]);
                    varsXplusDs.emplace_back(_secretData.XplusDs[ref]);

                    for (size_t offset : {0,1}) {
                        std::vector<FieldT> bits;
                        for (size_t bitIdx = 0; bitIdx < _secretData.bitDecompositions[2*idx+offset].size(); bitIdx++) {
                            bits.emplace_back(_secretData.bitDecompositions[2*idx+offset][bitIdx]);
                        }
                        assert(_secretData.bitDecompositions[2*idx+offset].size() == _publicData.log2Ds[idx]);
                        bitDecomposition.emplace_back(bits);
                    }
                }
            } else {
                varsXs.assign(sizeIndices, FieldT(0));
                varsXplusCs.assign(sizeIndices, FieldT(0));
                varsXplusDs.assign(sizeIndices, FieldT(0));

                for (size_t idx = 0; idx < sizeIndices; idx++) {
                    for (size_t offset : {0,1}) {
                        std::vector<FieldT> bits(_publicData.log2Ds[idx]);
                        bitDecomposition.emplace_back(bits);
                    }
                }
            }

            addIntrablocBoundingConstraintLarge(varsXs, varsXplusCs, varsXplusDs, Cs, ds, log2Ds, bitDecomposition, blocOneIdx, blocZeroIdx);
    }

    // All helper NP Expression functions
    void NP_RSARound3_keyGen(size_t blocZeroIdx, size_t blocOneIdx) {

            // Pick Evaluation Points Using Fiat-Shamir Transform
            std::vector<FieldT> evaluationPoints = fsPickEvaluationPoints();

            // Linear constraint Ax = b according to the Ligero paper, i.e. transformation according to the builder
            size_t blocCol = 0; // This variable keeps track of the current index for the misc bloc, in which we will mix a number of outputs

            transformation_constraint<FieldT> *transformationConstraint = new transformation_constraint<FieldT>;
            addKeyGenPolyEval(transformationConstraint, evaluationPoints, blocZeroIdx, blocOneIdx, blocCol);
            // padBloc(transformationConstraint, blocCol, blocZeroIdx);
            size_t keyGenEvalBlocIdx = _builder.add(transformationConstraint);

            // Then set this misc block equal to 0
            _builder.ensure_equality(keyGenEvalBlocIdx, blocZeroIdx);

    }

void NP_Equate_Variables(size_t blocZeroIdx, size_t blocOneIdx) {
        size_t blocCol = 0; // This variable keeps track of the current index for the misc bloc, in which we will mix a number of outputs

        std::vector<size_t> blocs_left;
        std::vector<size_t> cols_left;
        std::vector<size_t> blocs_right;
        std::vector<size_t> cols_right;

        // Populate with coordinates for the variables to be equated

        // Tying x_sharesPS and x_sharesPS_2
        blocs_left.insert(blocs_left.end(), this->_xsharesPSidxs.blocs.begin(), this->_xsharesPSidxs.blocs.end());
        cols_left.insert(cols_left.end(), this->_xsharesPSidxs.cols.begin(), this->_xsharesPSidxs.cols.end());
        blocs_right.insert(blocs_right.end(), this->_xsharesPSidxs_2.blocs.begin(), this->_xsharesPSidxs_2.blocs.end());
        cols_right.insert(cols_right.end(), this->_xsharesPSidxs_2.cols.begin(), this->_xsharesPSidxs_2.cols.end());

        // Tying x_sharesPS and x_sharesPS_3
        blocs_left.insert(blocs_left.end(), this->_xsharesPSidxs.blocs.begin(), this->_xsharesPSidxs.blocs.end());
        cols_left.insert(cols_left.end(), this->_xsharesPSidxs.cols.begin(), this->_xsharesPSidxs.cols.end());
        blocs_right.insert(blocs_right.end(), this->_xsharesPSidxs_3.blocs.begin(), this->_xsharesPSidxs_3.blocs.end());
        cols_right.insert(cols_right.end(), this->_xsharesPSidxs_3.cols.begin(), this->_xsharesPSidxs_3.cols.end());

        // Tying x_sharesPS and x_sharesPS_4
        blocs_left.insert(blocs_left.end(), this->_xsharesPSidxs.blocs.begin(), this->_xsharesPSidxs.blocs.end());
        cols_left.insert(cols_left.end(), this->_xsharesPSidxs.cols.begin(), this->_xsharesPSidxs.cols.end());
        blocs_right.insert(blocs_right.end(), this->_xsharesPSidxs_4.blocs.begin(), this->_xsharesPSidxs_4.blocs.end());
        cols_right.insert(cols_right.end(), this->_xsharesPSidxs_4.cols.begin(), this->_xsharesPSidxs_4.cols.end());

        // Tying y_sharesPS and y_sharesPS_2
        blocs_left.insert(blocs_left.end(), this->_ysharesPSidxs.blocs.begin(), this->_ysharesPSidxs.blocs.end());
        cols_left.insert(cols_left.end(), this->_ysharesPSidxs.cols.begin(), this->_ysharesPSidxs.cols.end());
        blocs_right.insert(blocs_right.end(), this->_ysharesPSidxs_2.blocs.begin(), this->_ysharesPSidxs_2.blocs.end());
        cols_right.insert(cols_right.end(), this->_ysharesPSidxs_2.cols.begin(), this->_ysharesPSidxs_2.cols.end());

        // Tying y_sharesPS and y_sharesPS_3
        blocs_left.insert(blocs_left.end(), this->_ysharesPSidxs.blocs.begin(), this->_ysharesPSidxs.blocs.end());
        cols_left.insert(cols_left.end(), this->_ysharesPSidxs.cols.begin(), this->_ysharesPSidxs.cols.end());
        blocs_right.insert(blocs_right.end(), this->_ysharesPSidxs_3.blocs.begin(), this->_ysharesPSidxs_3.blocs.end());
        cols_right.insert(cols_right.end(), this->_ysharesPSidxs_3.cols.begin(), this->_ysharesPSidxs_3.cols.end());

        // Tying y_sharesPS and y_sharesPS_4
        blocs_left.insert(blocs_left.end(), this->_ysharesPSidxs.blocs.begin(), this->_ysharesPSidxs.blocs.end());
        cols_left.insert(cols_left.end(), this->_ysharesPSidxs.cols.begin(), this->_ysharesPSidxs.cols.end());
        blocs_right.insert(blocs_right.end(), this->_ysharesPSidxs_4.blocs.begin(), this->_ysharesPSidxs_4.blocs.end());
        cols_right.insert(cols_right.end(), this->_ysharesPSidxs_4.cols.begin(), this->_ysharesPSidxs_4.cols.end());

        // Make sure all these variables tie out
        _builder.ensure_equality(
            addIntrablocEquality(blocZeroIdx, blocs_left, cols_left, blocs_right, cols_right), 
            blocZeroIdx);
}

void NP_RSARound4_preSieving(size_t blocZeroIdx, size_t blocOneIdx) {
        // ci_1 = u x a + v, ci_2 = u x b + w + (Q/P) m'
        // Pick Evaluation Points Using Fiat-Shamir Transform
        std::vector<FieldT> evaluationPoints = fsPickEvaluationPoints();

        // Linear constraint Ax = b according to the Ligero paper, i.e. transformation according to the builder
        size_t blocCol = 0; // This variable keeps track of the current index for the misc bloc, in which we will mix a number of outputs

        transformation_constraint<FieldT> *transformationConstraint = new transformation_constraint<FieldT>;

        // First Constraint: ci_1 = u * a + v
        //                      a = x * b + y
        addEquation1( 
            transformationConstraint,
            evaluationPoints, 
            blocZeroIdx, 
            blocOneIdx, 
            _publicData.ci_1, 
            _publicData.A,
            _secretData.ux,
            _secretData.vx,
            _secretData.q_r4_1,
            blocCol);

        // Second Constraint: ci_2 = u * b + w + (Q/P) m'
        //                       a = x * b + y + z
        addEquation2( 
            transformationConstraint, 
            evaluationPoints, 
            blocZeroIdx, 
            blocOneIdx, 
            _publicData.ci_2,
            _publicData.b,
            _secretData.ux,
            _secretData.wx,
            _secretData.xi,
            _secretData.q_r4_2,
            blocCol);

        // Padding the rest of the bloc
        // padBloc(transformationConstraint, blocCol, blocZeroIdx);
        size_t encryptionEvalBlocIdx = _builder.add(transformationConstraint);

        // Then set this misc block equal to 0
        _builder.ensure_equality(encryptionEvalBlocIdx, blocZeroIdx);

        if (_publicData.modulusIdx == 0) {
            // Bounding vi and wi
            if (_secretData.isAvailable) {
                addIntrablocBoundingConstraint(embed(_secretData.vxP,0,_secretData.vxP.size()), std::vector<FieldT>(_secretData.vxP.size(), FieldT(_publicData.sigma*10)), blocOneIdx, blocZeroIdx);
                addIntrablocBoundingConstraint(embed(_secretData.wxP,0,_secretData.wxP.size()), std::vector<FieldT>(_secretData.wxP.size(), FieldT(_publicData.sigma*10)), blocOneIdx, blocZeroIdx);
            } else {
                std::vector<FieldT> placeholder(30, FieldT(0));
                addIntrablocBoundingConstraint(placeholder, std::vector<FieldT>(placeholder.size(), FieldT(_publicData.sigma*10)), blocOneIdx, blocZeroIdx);
                addIntrablocBoundingConstraint(placeholder, std::vector<FieldT>(placeholder.size(), FieldT(_publicData.sigma*10)), blocOneIdx, blocZeroIdx);
            }
        }
    }

void NP_RSARound5_preSieving(size_t blocZeroIdx, size_t blocOneIdx) {
            // in this step:

            // encryption of zi (r5_ui, r5_vi, r5_wi)
            // linear constraint for ci_1_prime[1] and ci_1_prime[2]

            // Pick Evaluation Points Using Fiat-Shamir Transform
            std::vector<FieldT> evaluationPoints = fsPickEvaluationPoints();

            // Linear constraint Ax = b according to the Ligero paper, i.e. transformation according to the builder
            size_t blocCol = 0; // This variable keeps track of the current index for the misc bloc, in which we will mix a number of outputs

            transformation_constraint<FieldT> *transformationConstraint = new transformation_constraint<FieldT>;

            // First Constraint: ci_1 = yi*c_1 -u*a - v
            //                      a = x*b    -y*c - z

            addEquation3( 
                transformationConstraint,
                evaluationPoints,
                blocZeroIdx,
                blocOneIdx,
                _publicData.ci_1_prime,
                _publicData.c_1,
                _publicData.A,
                _secretData.yi,
                _secretData.uz,
                _secretData.vz,
                _secretData.q_r5_1,
                blocCol);

            // Second Constraint: ci_2 = yi*c_2 + zp - uz*b - wz - (Q/P) m'
            //                       a = x*b    + z  - y*c  - t  -  m

            addEquation4( 
                transformationConstraint,
                evaluationPoints,
                blocZeroIdx,
                blocOneIdx,
                _publicData.ci_2_prime,
                _publicData.c_2,
                _publicData.b,
                _secretData.yi,
                _secretData.uz,
                _secretData.zp, 
                _secretData.wz,
                _secretData.zi, 
                _secretData.q_r5_2, // qci_2
                blocCol);

            // Padding the rest of the bloc
            // padBloc(transformationConstraint, blocCol, blocZeroIdx);
            size_t encryptionEvalBlocIdx = _builder.add(transformationConstraint);

            // Then set this misc block equal to 0
            _builder.ensure_equality(encryptionEvalBlocIdx, blocZeroIdx);

            if (_publicData.modulusIdx == 0) {
                // Bounding vi and wi
                if (_secretData.isAvailable) {
                    addIntrablocBoundingConstraint(embed(_secretData.vzP,0,_secretData.vzP.size()), std::vector<FieldT>(_secretData.vzP.size(), FieldT(_publicData.sigma*10)), blocOneIdx, blocZeroIdx);
                    addIntrablocBoundingConstraint(embed(_secretData.wzP,0,_secretData.wzP.size()), std::vector<FieldT>(_secretData.wzP.size(), FieldT(_publicData.sigma*10)), blocOneIdx, blocZeroIdx);
                } else {
                    std::vector<FieldT> placeholder(30, FieldT(0));
                    addIntrablocBoundingConstraint(placeholder, std::vector<FieldT>(placeholder.size(), FieldT(_publicData.sigma*10)), blocOneIdx, blocZeroIdx);
                    addIntrablocBoundingConstraint(placeholder, std::vector<FieldT>(placeholder.size(), FieldT(_publicData.sigma*10)), blocOneIdx, blocZeroIdx);
                }
            }
    }

void NP_RSARound6_partialDecrypt(size_t blocZeroIdx, size_t blocOneIdx) {
        // in this step:
        // partial decryption

        // Pick Evaluation Points Using Fiat-Shamir Transform
        std::vector<FieldT> evaluationPoints = fsPickEvaluationPoints();

        // Linear constraint Ax = b according to the Ligero paper, i.e. transformation according to the builder
        size_t blocCol = 0; // This variable keeps track of the current index for the misc bloc, in which we will mix a number of outputs

        transformation_constraint<FieldT> *transformationConstraint = new transformation_constraint<FieldT>;

        // First Constraint:   di = c'_2 - c'_1*s + r;
        //                      a = c    - x*b    + y 
        addEquation5( 
            transformationConstraint,
            evaluationPoints,
            blocZeroIdx,
            blocOneIdx,
            _publicData.di,
            _publicData.c_1_prime,
            _publicData.c_2_prime,
            _secretData.siP,
            _secretData.rNoise,
            _secretData.q_r6,
            blocCol);

        // Padding the rest of the bloc
        // padBloc(transformationConstraint, blocCol, blocZeroIdx);
        size_t encryptionEvalBlocIdx = _builder.add(transformationConstraint);

        // Then set this misc block equal to 0
        _builder.ensure_equality(encryptionEvalBlocIdx, blocZeroIdx);
    }

void checkEvaluation(transformation_constraint<FieldT> *tc, std::vector<size_t> coefsBlocIdxs, std::vector<size_t> points, std::vector<size_t>& blocs, std::vector<size_t>& cols, std::vector<uint64_t>& qdivp, size_t& blocCol, bool scaling) {
        size_t i = 0;

        for (size_t eval : points) {
            std::vector<FieldT> powersOfEval(_publicData.degree);
            FieldT tmp(1);

            for (size_t powIdx = 0; powIdx < _publicData.degree; powIdx++) {
                powersOfEval[powIdx] = tmp;
                tmp = tmp * this->_phis[eval];
            }

            // Adding f(x) = y constraint
            assert (blocCol < this->_blockSize);

            // debug: check that this evaluation works
            #ifndef NDEBUG
            if (_secretData.isAvailable) {
                nfl::poly_p<uint64_t, 65536, 21> xi;
                for (size_t modIdx = 0; modIdx < 21; modIdx++) {
                    for (size_t powIdx = 0; powIdx < _publicData.degree; powIdx++) {
                        xi.poly_obj().data()[powIdx+modIdx*_publicData.degree] = this->_secretData.xi[powIdx+modIdx*_publicData.degree];
                    }
                }

                xi.ntt_pow_phi();

                FieldT evalManual(0);

                for (size_t powIdx = 0; powIdx < _publicData.degree; powIdx++) {
                    evalManual = evalManual + powersOfEval[powIdx]*FieldT(this->_secretData.xi[powIdx+_publicData.modulusIdx*_publicData.degree]);
                }

                // std::cout << "evalManual" << evalManual.getValue() << std::endl;
                // std::cout << "xi:" << xi.poly_obj().data()[_publicData.indicesGCD[i]] << std::endl;
                // std::cout << "x_shares:" << this->_secretData.x_sharesGCD[i] << std::endl;
                // std::cout << "q_over_p:" << this->_qdivpPS[i] << std::endl;
                // std::cout << "modulus:" << FieldT::getModulus() << std::endl;

                // std::cout << evalManual.getValue() << std::endl;
                // std::cout << static_cast<uint64_t>(static_cast<__uint128_t>(this->_secretData.x_sharesGCD[i])*static_cast<__uint128_t>(this->_qdivpGCD[i]) % static_cast<__uint128_t>(FieldT::getModulus())) << std::endl;

                // std::cout << (evalManual-(FieldT(this->_secretData.x_sharesGCD[i*21])*FieldT(this->_qdivpGCD[i]))).getValue() << std::endl;
            }
            #endif

            addEvalPt(tc, blocCol, FieldT(1), coefsBlocIdxs, powersOfEval);
            if (scaling) connectValuations(tc, blocs, cols, qdivp, blocCol, i);
            else connectValuations(tc, blocs, cols, blocCol, i);

            i++;blocCol++;
        }
}

void connectValuations(transformation_constraint<FieldT> *tc, std::vector<size_t>& blocs, std::vector<size_t>& cols, size_t& blocCol, size_t idx) {
        FieldT minusOne(FieldT::getModulus()-1);

        tc->block.emplace_back(blocs[idx]);
        tc->source_position.emplace_back(cols[idx]);
        tc->scalar.emplace_back(minusOne);
        tc->target_position.emplace_back(blocCol);
}

void connectValuations(transformation_constraint<FieldT> *tc, std::vector<size_t>& blocs, std::vector<size_t>& cols, std::vector<uint64_t>& qdivp, size_t& blocCol, size_t idx) {
        FieldT minusOne(FieldT::getModulus()-1);

        tc->block.emplace_back(blocs[idx]);
        tc->source_position.emplace_back(cols[idx]);
        tc->scalar.emplace_back(FieldT(qdivp[idx]) * minusOne);
        tc->target_position.emplace_back(blocCol);
}

void NP_Connecting_Proofs(size_t blocZeroIdx, size_t blocOneIdx) {
        size_t blocCol = 0;
        transformation_constraint<FieldT> *tc = new transformation_constraint<FieldT>;

        // build evaluation points from PS, CAN and GCD
        checkEvaluation(tc, this->_xprimeIds,this->_publicData.indicesPS, this->_xsharesPSidxs.blocs, this->_xsharesPSidxs.cols, this->_qdivpPS, blocCol, true);
        checkEvaluation(tc, this->_xprimeIds,this->_publicData.indicesCAN, this->_xsharesCANidxs.blocs, this->_xsharesCANidxs.cols, this->_qdivpCAN, blocCol, true);
        checkEvaluation(tc, this->_xprimeIds,this->_publicData.indicesGCD, this->_xsharesGCDidxs.blocs, this->_xsharesGCDidxs.cols, this->_qdivpGCD, blocCol, true);

        checkEvaluation(tc, this->_yprimeIds,this->_publicData.indicesPS, this->_ysharesPSidxs.blocs, this->_ysharesPSidxs.cols, this->_qdivpPS, blocCol, false);
        checkEvaluation(tc, this->_yprimeIds,this->_publicData.indicesCAN, this->_ysharesCANidxs.blocs, this->_ysharesCANidxs.cols, this->_qdivpCAN, blocCol, false);
        checkEvaluation(tc, this->_yprimeIds,this->_publicData.indicesGCD, this->_ysharesGCDidxs.blocs, this->_ysharesGCDidxs.cols, this->_qdivpGCD, blocCol, false);

        checkEvaluation(tc, this->_zprimeIds,this->_publicData.indicesCAN, this->_zsharesCANidxs.blocs, this->_zsharesCANidxs.cols, this->_qdivpCAN, blocCol, true);

        // Padding the rest of the bloc
        // padBloc(tc, blocCol, blocZeroIdx);
        size_t connectingProofsRowIdx = _builder.add(tc);

        // Then set this block equal to 0
        _builder.ensure_equality(connectingProofsRowIdx, blocZeroIdx);
}

void NP_RSARound7_candidateGeneration(size_t blocZeroIdx, size_t blocOneIdx) {

        // in this step:
        // partial decryption

        // Linear constraint Ax = b according to the Ligero paper, i.e. transformation according to the builder
        size_t blocCol = 0; // This variable keeps track of the current index for the misc bloc, in which we will mix a number of outputs

        transformation_constraint<FieldT> *tc = new transformation_constraint<FieldT>;

        auto [x_sharesPSidxs, x_sharesCANidxs] = addEquation0(
            tc,
            blocZeroIdx,
            blocOneIdx,
            _secretData.x_sharesPS,
            _secretData.x_sharesCAN,
            _secretData.q_p_prod_r7,
            _secretData.q_p_r7,
            _publicData.ax_shares,
            _publicData.coefsCAN,
            _publicData.prodcans,
            _publicData.cans,
            blocCol
        );

        this->_xsharesPSidxs = x_sharesPSidxs;
        this->_xsharesCANidxs = x_sharesCANidxs;

        auto [y_sharesPSidxs, y_sharesCANidxs] = addEquation0(
            tc,
            blocZeroIdx,
            blocOneIdx,
            _secretData.y_sharesPS,
            _secretData.y_sharesCAN,
            _secretData.q_q_prod_r7,
            _secretData.q_q_r7,
            _publicData.by_shares,
            _publicData.coefsCAN,
            _publicData.prodcans,
            _publicData.cans,
            blocCol
        );

        this->_ysharesPSidxs = y_sharesPSidxs;
        this->_ysharesCANidxs = y_sharesCANidxs;

        // Padding the rest of the bloc
        // padBloc(tc, blocCol, blocZeroIdx);
        size_t canGenEvalBlocIdx = _builder.add(tc);

        // Then set this misc block equal to 0
        _builder.ensure_equality(canGenEvalBlocIdx, blocZeroIdx);
    }

    void NP_RSARound8_beaversTriples(size_t blocZeroIdx, size_t blocOneIdx) {
            // Linear constraint Ax = b according to the Ligero paper, i.e. transformation according to the builder
            size_t blocCol = 0; // This variable keeps track of the current index for the misc bloc, in which we will mix a number of outputs

            transformation_constraint<FieldT> *tc = new transformation_constraint<FieldT>;

            addEquation00(
                tc,
                blocZeroIdx, 
                blocOneIdx,
                _secretData.y_sharesPS,
                _secretData.q_q_prod_r7,
                _publicData.ax,
                _publicData.coefsCAN,
                _publicData.prodcans,

                _secretData.x_sharesPS,
                _secretData.q_p_prod_r7,
                _publicData.by,
                _publicData.prodcans,
                
                _secretData.z_sharesCAN,
                _secretData.q_r8,
                _publicData.axby,
                _publicData.cans,
                blocCol
            );

            // Padding the rest of the bloc
            // padBloc(tc, blocCol, blocZeroIdx);
            size_t canGenEvalBlocIdx = _builder.add(tc);

            // Then set this misc block equal to 0
            _builder.ensure_equality(canGenEvalBlocIdx, blocZeroIdx);
    }

    void NP_RSARound10_jacobiTest(size_t blocZeroIdx, size_t blocOneIdx) {}

    void NP_RSARound11_12_jacobiGCD(size_t blocZeroIdx, size_t blocOneIdx) {
            // Linear constraint Ax = b according to the Ligero paper, i.e. transformation according to the builder
            size_t blocCol = 0; // This variable keeps track of the current index for the misc bloc, in which we will mix a number of outputs

            transformation_constraint<FieldT> *tc = new transformation_constraint<FieldT>;

            addEquation11(
                tc,
                blocZeroIdx,
                blocOneIdx,

                //0 [11]
                _secretData.y_sharesGCD,
                _publicData.by_shares_GCD,

                //1
                _secretData.x_sharesPS,
                _secretData.q_p_prod_r11,
                _publicData.coefsGCD,
                _publicData.prodgcds,

                //2
                _secretData.y_sharesPS,
                _secretData.q_q_prod_r11,
                _secretData.q_pq_r11,
                _publicData.coefsGCD,
                _publicData.prodgcds,
                _publicData.gcds,

                //3 [12]
                _secretData.x_sharesGCD,
                _secretData.r_CRTs,
                _secretData.q_r_r11,
                _publicData.ax_shares_GCD,
                _publicData.gcds,
                blocCol
            );

            addEquation12(
                tc,
                blocZeroIdx,
                blocOneIdx,

                //0 [11]
                _publicData.axbyGCD,

                //1
                _secretData.x_sharesPS,
                _secretData.q_p_prod_r11,
                _secretData.y_sharesPS,
                _secretData.q_q_prod_r11,
                _publicData.axGCD,
                _publicData.coefsGCD,
                _publicData.prodgcds,

                //3 [12]
                _secretData.r_CRTs,
                _secretData.z_sharesGCD,
                _secretData.q_r12,
                _publicData.byGCD,
                _publicData.gcds,
                blocCol
            );

            // Padding the rest of the bloc
            // padBloc(tc, blocCol, blocZeroIdx);
            size_t canGenEvalBlocIdx = _builder.add(tc);

            // Then set this misc block equal to 0
            _builder.ensure_equality(canGenEvalBlocIdx, blocZeroIdx);
    }

    // For evaluation purpose only: in practice this will be too onerous for use
    std::vector<size_t> addConvolutionConstraint(std::vector<size_t> sourceIdx, size_t evaluationDomainSize) {

        // var: nShares
        DBG("NFLlib Prime:" << FieldT::getModulus());
        DBG("NFLlib Root w:" << _publicData.roots[_publicData.modulusIdx*_publicData.degree]);

        // Add intrabloc constraints
        std::vector<size_t> convolutionTargetBlocsIdx;
        transformation_constraint<FieldT> *newTransformationConstraint = new transformation_constraint<FieldT>[evaluationDomainSize*evaluationDomainSize];
        size_t idx = 0;
        
        for (size_t strideBloc = 0; strideBloc < (evaluationDomainSize/this->_blockSize); strideBloc++) {
            for (size_t stride = 0; stride < _blockSize; stride++) {
                for (size_t step = 0; step < evaluationDomainSize; step++) {
                    // Source the adequate coefficients
                    newTransformationConstraint[idx].block.emplace_back(sourceIdx[step/this->_blockSize]);
                    newTransformationConstraint[idx].source_position.emplace_back(step % this->_blockSize);

                    size_t evaluationDomainIdx = strideBloc*this->_blockSize + stride;
                    FieldT evalPt(_publicData.roots[evaluationDomainIdx]);
                    FieldT eval = (evalPt^step).getValue();

                    newTransformationConstraint[idx].scalar.emplace_back(eval);
                    newTransformationConstraint[idx].target_position.emplace_back(stride);
                }
            }
            convolutionTargetBlocsIdx.emplace_back(_builder.add(&newTransformationConstraint[idx++]));
        }

        return convolutionTargetBlocsIdx;
    }

    size_t addPolyEvalConstraint(std::vector<size_t>& sourceIdx, size_t evaluationDomainSize, std::vector<size_t> positions, size_t zeroBlocRef) {

        // var: nShares
        DBG("NFLlib Prime:" << FieldT::getModulus());
        DBG("NFLlib Root w:" << _publicData.roots[0]);

        // Add intrabloc constraints
        transformation_constraint<FieldT> *newTransformationConstraint = new transformation_constraint<FieldT>;

        DBG("size of positions:" << positions.size());

        // Regardless of the actual position observed, we will always place it at the begining of the bloc
        size_t blocCol = 0;
        for (size_t idxEval : positions) {
            for (size_t step = 0; step < evaluationDomainSize; step++) {
                // Source the adequate coefficients
                newTransformationConstraint->block.emplace_back(sourceIdx[step/this->_blockSize]);
                newTransformationConstraint->source_position.emplace_back(step % this->_blockSize);

                FieldT evalPt(_publicData.roots[idxEval]);
                FieldT eval = (evalPt^step).getValue();

                newTransformationConstraint->scalar.emplace_back(eval);
                newTransformationConstraint->target_position.emplace_back(blocCol);
            }
            
            blocCol++;
        }

        // Pad the rest of the bloc
        for (size_t idxEval = positions.size(); idxEval < this->_blockSize; idxEval++) {
            newTransformationConstraint->block.emplace_back(zeroBlocRef);
            newTransformationConstraint->source_position.emplace_back(0);
            newTransformationConstraint->scalar.emplace_back(0);
            newTransformationConstraint->target_position.emplace_back(idxEval);
        }

        return _builder.add(newTransformationConstraint);
    }

    // Adding the evaluation of a polynomial on a specific point
    void addEvalPt(transformation_constraint<FieldT> *tc, size_t targetPosition, FieldT scalar, std::vector<size_t>& coefRefs, std::vector<FieldT>& powerEval) { 

        for (size_t step = 0; step < _publicData.degree; step++) {
            tc->block.emplace_back(coefRefs[step/this->_blockSize]);
            tc->source_position.emplace_back(step % this->_blockSize);
            tc->scalar.emplace_back(powerEval[step]*scalar);
            tc->target_position.emplace_back(targetPosition);
        }
    }

    // Here adding a simple scalar on a specific point
    void addEvalPt(transformation_constraint<FieldT> *tc, size_t targetPosition, FieldT scalar, size_t blocOneIdxRef) { 

        tc->block.emplace_back(blocOneIdxRef);
        tc->source_position.emplace_back(0);
        tc->scalar.emplace_back(scalar);
        tc->target_position.emplace_back(targetPosition);
    }

    // Simple multiplication
    // TGT = sum ( A*X + Y - C - B*Q )
    void addMults(transformation_constraint<FieldT> *tc, size_t sourceBlock, size_t targetPosition, std::vector<FieldT>& scalarsA, std::vector<size_t>& positionsX, size_t positionY, FieldT scalarC, FieldT scalarB, size_t positionQ, size_t blockOneIdx) { 

        // sum (x*a)
        for (size_t idx = 0; idx < scalarsA.size(); idx++) {
            tc->block.emplace_back(sourceBlock);
            tc->source_position.emplace_back(positionsX[idx]);
            tc->scalar.emplace_back(scalarsA[idx]);
            tc->target_position.emplace_back(targetPosition);
        }

        // + y
        tc->block.emplace_back(sourceBlock);
        tc->source_position.emplace_back(positionY);
        tc->scalar.emplace_back(FieldT(1));
        tc->target_position.emplace_back(targetPosition);

        // - c
        tc->block.emplace_back(blockOneIdx);
        tc->source_position.emplace_back(0);
        tc->scalar.emplace_back(FieldT(FieldT::getModulus()-1)*scalarC);
        tc->target_position.emplace_back(targetPosition);

        // - q*b
        tc->block.emplace_back(sourceBlock);
        tc->source_position.emplace_back(positionQ);
        tc->scalar.emplace_back(scalarB);
        tc->target_position.emplace_back(targetPosition);
    }

    FieldT naiveEval(std::vector<FieldT> coefs, std::vector<FieldT> powersOfEval) {
        FieldT returnValue(0);

        for (size_t idx = 0; idx < coefs.size(); idx++) {
            returnValue += coefs[idx]*powersOfEval[idx];
        }

        return returnValue;
    }

    void addKeyGenPolyEval(transformation_constraint<FieldT> *newTransformationConstraint, std::vector<FieldT>& evaluationPoints, size_t blocZeroIdx, size_t blockOneIdx, size_t& blocCol) {

        // We assume that the bloc size divides the poly evaluation size for simplicity
        size_t nbBlocs = (_publicData.degree/_blockSize);    

        // Source All Data
        // ======================================================================================================

        // Public Data
        // A and bi are stored in the evaluation domain
        std::vector<FieldT> a = embed(_publicData.A, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        std::vector<FieldT> bi = embed(_publicData.bi, _publicData.modulusIdx * _publicData.degree, _publicData.degree);

        std::vector<FieldT> ei;
        std::vector<FieldT> si;
        std::vector<FieldT> qi;

        // Secret Data: either source or initialize if not available
        if (_secretData.isAvailable) {
            ei = embed(_secretData.eiP, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            si = embed(_secretData.siP, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            qi = embed(_secretData.q_r3, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        } else {
            ei.assign(_publicData.degree, FieldT(0));
            si.assign(_publicData.degree, FieldT(0));
            qi.assign(_publicData.degree, FieldT(0));
        }

        // Building the constraint system
        // ======================================================================================================

        // Bounding Variables

        // Bounding Constraint
        // Rephrase the problem using the nearest 2^n
        // Here we look to bound by 10*sigma, so we need D larger than 2* (10*sigma) 
        // to phrase the problem as absolute numbers
        auto [d, C, D, log2D] = rephraseAs2Nct(_publicData.sigma*10);

        vector<size_t> siCoefsBlocIds, eiCoefsBlocIds;

        if (_publicData.modulusIdx ==0) {
            // Bounding variables
            siCoefsBlocIds = addBoundingConstraint(_publicData.degree, si, FieldT(C), FieldT(d), static_cast<uint64_t>(log2D), blockOneIdx);
            eiCoefsBlocIds = addBoundingConstraint(_publicData.degree, ei, FieldT(C), FieldT(d), static_cast<uint64_t>(log2D), blockOneIdx);
            assert(siCoefsBlocIds.size() == eiCoefsBlocIds.size());
        } else {
            siCoefsBlocIds = addSplitVariable(si);
            eiCoefsBlocIds = addSplitVariable(ei);
        }

        // Splitting variables
        vector<size_t> qiCoefsBlocIds = addSplitVariable(qi);

        DBG("Number of evaluation points:" << evaluationPoints.size());

        for (FieldT eval : evaluationPoints) {
            std::vector<FieldT> powersOfEval(_publicData.degree);

            for (size_t powIdx = 0; powIdx < _publicData.degree; powIdx++) {
                powersOfEval[powIdx] = eval^powIdx;
            }

            DBG("Adding a[j]*si[j] + e[j] - b[j] == (x^n + 1) * q[j] constraint:");
            assert (blocCol < this->_blockSize);

            addEvalPt(newTransformationConstraint, blocCol, naiveEval(a, powersOfEval), siCoefsBlocIds, powersOfEval); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(1), eiCoefsBlocIds, powersOfEval); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1)*naiveEval(bi, powersOfEval), blockOneIdx); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1)*(eval*powersOfEval.back()+1), qiCoefsBlocIds, powersOfEval); 

            blocCol++;
        }
    }

    void padBloc(transformation_constraint<FieldT> *newTransformationConstraint, size_t &blocCol, size_t blocZeroIdx) {

        // Pad the rest of the bloc
        for (; blocCol < this->_blockSize; blocCol++) {
            newTransformationConstraint->block.emplace_back(blocZeroIdx);
            newTransformationConstraint->source_position.emplace_back(0);
            newTransformationConstraint->scalar.emplace_back(0);
            newTransformationConstraint->target_position.emplace_back(blocCol);
        }
    }

    auto addEquation0(
        transformation_constraint<FieldT> *tc, 
        size_t blocZeroIdx,
        size_t blocOneIdx,
        std::vector<uint64_t>& x,
        std::vector<uint64_t>& y,
        std::vector<uint64_t>& z,
        std::vector<uint64_t>& t,
        std::vector<uint64_t>& a,
        std::vector<uint64_t>& b,
        std::vector<uint64_t>& c,
        std::vector<uint64_t>& d,
        size_t& blocCol
    ) {
        // populating the misc bloc

        size_t& p = _publicData.modulusIdx;

        // per prime
        std::vector<FieldT> miscBloc(this->_blockSize, FieldT(0));
        size_t miscIdx=0;

        auto addEmptyVar = [&]() -> size_t {
            miscBloc[miscIdx++] = FieldT(0);
            return (miscIdx-1);
        };

        auto addEmptyVect = [&](std::vector<size_t> idxs) -> size_t {
            for (size_t idx : idxs) {miscBloc[miscIdx++] = FieldT(0);}
            return (miscIdx - idxs.size());
        };

        auto addVar = [&](uint64_t a) -> size_t {
            miscBloc[miscIdx++] = a;
            return (miscIdx-1);
        };

        auto addVect = [&](std::vector<uint64_t>& a, std::vector<size_t> idxs) -> size_t {
            for (size_t idx : idxs) {miscBloc[miscIdx++] = a[idx];}
            return (miscIdx - idxs.size());
        };

        // Declaring all variables
        size_t s_x_0;
        size_t s_x_1;
        size_t s_x_2;
        size_t s_x_3;
        size_t s_x_4;
        size_t s_x_5;

        size_t s_y_idx;
        size_t s_z_idx;
        size_t s_t_idx;

        std::pair<intraBlocReferences, intraBlocReferences> ret;

        // Vectors
        std::vector<size_t> v;
        for (size_t idx = 0; idx < alphaCanLength; idx++) {v.push_back((idx*NbPrimesQ)+p);}

        if (_secretData.isAvailable) {

            // Single Variables
            s_x_0   = addVar(x[              p]);
            s_x_1   = addVar(x[  NbPrimesQ + p]);
            s_x_2   = addVar(x[2*NbPrimesQ + p]);
            s_x_3   = addVar(x[3*NbPrimesQ + p]);
            s_x_4   = addVar(x[4*NbPrimesQ + p]);
            s_x_5   = addVar(x[5*NbPrimesQ + p]);

            s_y_idx    = addVect(y, v);
            s_z_idx    = addVect(z, v);
            s_t_idx    = addVect(t, v);

        } else {
            // Single Variables
            s_x_0   = addEmptyVar();
            s_x_1   = addEmptyVar();
            s_x_2   = addEmptyVar();
            s_x_3   = addEmptyVar();
            s_x_4   = addEmptyVar();
            s_x_5   = addEmptyVar();

            s_y_idx    = addEmptyVect(v);
            s_z_idx    = addEmptyVect(v);
            s_t_idx    = addEmptyVect(v);
        }

        // Shared book-keeping for connecting proofs
        ret.first.cols.insert(ret.first.cols.end(), {s_x_0,s_x_1,s_x_2,s_x_3,s_x_4,s_x_5});
        ret.second.cols.insert(ret.second.cols.end(), {s_y_idx,s_y_idx+1,s_y_idx+2,s_y_idx+3,s_y_idx+4,s_y_idx+5});

        // Populating inputs 
        size_t miscBlocIdx =_builder.add(new variable_constraint<FieldT>());
        _variables.add(miscBlocIdx,miscBloc);

        // Shared book-keeping for connecting proofs
        ret.first.blocs.assign(6, miscBlocIdx);
        ret.second.blocs.assign(6, miscBlocIdx);

        FieldT minusOne(FieldT::getModulus()-1);
        // Linear Combination
        for (size_t idx = 0; idx < alphaCanLength; idx++) {

            // y[i+p] + a[i+p]
            tc->block.emplace_back(blocOneIdx);
            tc->source_position.emplace_back(0);
            tc->scalar.emplace_back(FieldT(a[idx*NbPrimesQ +p]));
            tc->target_position.emplace_back(blocCol);

            tc->block.emplace_back(miscBlocIdx);
            tc->source_position.emplace_back(s_y_idx+idx);
            tc->scalar.emplace_back(FieldT(1));
            tc->target_position.emplace_back(blocCol);

            // 
            tc->block.emplace_back(miscBlocIdx);
            tc->source_position.emplace_back(s_x_0);
            tc->scalar.emplace_back(minusOne * FieldT(b[7*(idx*NbPrimesQ)+0*NbPrimesQ+p]));
            tc->target_position.emplace_back(blocCol);

            tc->block.emplace_back(miscBlocIdx);
            tc->source_position.emplace_back(s_x_1);
            tc->scalar.emplace_back(minusOne * FieldT(b[7*(idx*NbPrimesQ)+1*NbPrimesQ+p]));
            tc->target_position.emplace_back(blocCol);

            tc->block.emplace_back(miscBlocIdx);
            tc->source_position.emplace_back(s_x_2);
            tc->scalar.emplace_back(minusOne * FieldT(b[7*(idx*NbPrimesQ)+2*NbPrimesQ+p]));
            tc->target_position.emplace_back(blocCol);

            tc->block.emplace_back(miscBlocIdx);
            tc->source_position.emplace_back(s_x_3);
            tc->scalar.emplace_back(minusOne * FieldT(b[7*(idx*NbPrimesQ)+3*NbPrimesQ+p]));
            tc->target_position.emplace_back(blocCol);

            tc->block.emplace_back(miscBlocIdx);
            tc->source_position.emplace_back(s_x_4);
            tc->scalar.emplace_back(minusOne * FieldT(b[7*(idx*NbPrimesQ)+4*NbPrimesQ+p]));
            tc->target_position.emplace_back(blocCol);

            tc->block.emplace_back(miscBlocIdx);
            tc->source_position.emplace_back(s_x_5);
            tc->scalar.emplace_back(minusOne * FieldT(b[7*(idx*NbPrimesQ)+5*NbPrimesQ+p]));
            tc->target_position.emplace_back(blocCol);

            tc->block.emplace_back(blocOneIdx);
            tc->source_position.emplace_back(0);
            tc->scalar.emplace_back(minusOne * FieldT(_publicData.special ? FieldT(3) : FieldT(0))
                            * FieldT(b[7*(idx*NbPrimesQ)+6*NbPrimesQ+p]));
            tc->target_position.emplace_back(blocCol);

            tc->block.emplace_back(miscBlocIdx);
            tc->source_position.emplace_back(s_z_idx + idx);
            tc->scalar.emplace_back(FieldT(c[(idx*NbPrimesQ) +p]));
            tc->target_position.emplace_back(blocCol);

            tc->block.emplace_back(miscBlocIdx);
            tc->source_position.emplace_back(s_t_idx + idx);
            tc->scalar.emplace_back(minusOne * FieldT(d[(idx*NbPrimesQ) +p]));
            tc->target_position.emplace_back(blocCol);

            blocCol++;
        }

        return ret;
    }

    void addEquation00(
        transformation_constraint<FieldT> *tc, 
        size_t blocZeroIdx,
        size_t blocOneIdx,
        std::vector<uint64_t>& x,
        std::vector<uint64_t>& y,
        std::vector<uint64_t>& a,
        std::vector<uint64_t>& b,
        std::vector<uint64_t>& c,

        std::vector<uint64_t>& x_2,
        std::vector<uint64_t>& y_2,
        std::vector<uint64_t>& a_2,
        std::vector<uint64_t>& c_2,

        std::vector<uint64_t>& x_3,
        std::vector<uint64_t>& y_3,
        std::vector<uint64_t>& a_3,
        std::vector<uint64_t>& b_3,
        size_t& blocCol
    ) {
        // populating the misc bloc

        size_t& p = _publicData.modulusIdx;

        // per prime
        std::vector<FieldT> miscBloc(this->_blockSize, FieldT(0));
        size_t miscIdx=0;

        auto addEmptyVar = [&]() -> size_t {
            miscBloc[miscIdx++] = FieldT(0);
            return (miscIdx-1);
        };

        auto addEmptyVect = [&](std::vector<size_t> idxs) -> size_t {
            for (size_t idx : idxs) {miscBloc[miscIdx++] = FieldT(0);}
            return (miscIdx - idxs.size());
        };

        auto addVar = [&](uint64_t a) -> size_t {
            miscBloc[miscIdx++] = a;
            return (miscIdx-1);
        };

        auto addVect = [&](std::vector<uint64_t>& a, std::vector<size_t> idxs) -> size_t {
            for (size_t idx : idxs) {miscBloc[miscIdx++] = a[idx];}
            return (miscIdx - idxs.size());
        };

        constexpr int NbPrimesQ         = 21;
        constexpr int alphaCanLength    = 6;
        constexpr int alphaPSLength     = 6;

        // Declaring all variables
        size_t s_x_0;
        size_t s_x_1;
        size_t s_x_2;
        size_t s_x_3;
        size_t s_x_4;
        size_t s_x_5;

        size_t s_y_idx;

        size_t s_x_2_0;
        size_t s_x_2_1;
        size_t s_x_2_2;
        size_t s_x_2_3;
        size_t s_x_2_4;
        size_t s_x_2_5;

        size_t s_y_2_idx;

        size_t s_x_3_idx;
        size_t s_y_3_idx;

        // Vectors
        std::vector<size_t> v;
        for (size_t idx = 0; idx < alphaCanLength; idx++) {v.push_back(idx*NbPrimesQ+p);}

        if (_secretData.isAvailable) {

            // Single Variables
            s_x_0   = addVar(x[              p]);
            s_x_1   = addVar(x[  NbPrimesQ + p]);
            s_x_2   = addVar(x[2*NbPrimesQ + p]);
            s_x_3   = addVar(x[3*NbPrimesQ + p]);
            s_x_4   = addVar(x[4*NbPrimesQ + p]);
            s_x_5   = addVar(x[5*NbPrimesQ + p]);

            s_x_2_0   = addVar(x_2[              p]);
            s_x_2_1   = addVar(x_2[  NbPrimesQ + p]);
            s_x_2_2   = addVar(x_2[2*NbPrimesQ + p]);
            s_x_2_3   = addVar(x_2[3*NbPrimesQ + p]);
            s_x_2_4   = addVar(x_2[4*NbPrimesQ + p]);
            s_x_2_5   = addVar(x_2[5*NbPrimesQ + p]);

            s_y_idx      = addVect(y, v);
            s_y_2_idx    = addVect(y_2, v);
            
            s_x_3_idx    = addVect(x_3, v);
            s_y_3_idx    = addVect(y_3, v);

        } else {
            // Single Variables
            s_x_0   = addEmptyVar();
            s_x_1   = addEmptyVar();
            s_x_2   = addEmptyVar();
            s_x_3   = addEmptyVar();
            s_x_4   = addEmptyVar();
            s_x_5   = addEmptyVar();

            s_x_2_0   = addEmptyVar();
            s_x_2_1   = addEmptyVar();
            s_x_2_2   = addEmptyVar();
            s_x_2_3   = addEmptyVar();
            s_x_2_4   = addEmptyVar();
            s_x_2_5   = addEmptyVar();

            s_y_idx      = addEmptyVect(v);
            s_y_2_idx    = addEmptyVect(v);

            s_x_3_idx    = addEmptyVect(v);
            s_y_3_idx    = addEmptyVect(v);
        }

        // Book-keep location of x shares
        _xsharesPSidxs_2.cols.assign({s_x_2_0, s_x_2_1, s_x_2_2, s_x_2_3, s_x_2_4, s_x_2_5});
        _ysharesPSidxs_2.cols.assign({s_x_0, s_x_1, s_x_2, s_x_3, s_x_4, s_x_5});

        // Shared book-keeping for connecting proofs
        _zsharesCANidxs.cols.insert(_zsharesCANidxs.cols.end(), {s_x_3_idx,s_x_3_idx+1,s_x_3_idx+2,s_x_3_idx+3,s_x_3_idx+4,s_x_3_idx+5});

        // Populating inputs 
        size_t miscBlocIdx =_builder.add(new variable_constraint<FieldT>());
        _variables.add(miscBlocIdx,miscBloc);

        _zsharesCANidxs.blocs.assign(alphaCanLength, miscBlocIdx);
        _xsharesPSidxs_2.blocs.assign(alphaPSLength, miscBlocIdx);
        _ysharesPSidxs_2.blocs.assign(alphaPSLength, miscBlocIdx);

        FieldT minusOne(FieldT::getModulus()-1);
        // Linear Combination
        for (size_t idx = 0; idx < alphaCanLength; idx++) {

            {
                FieldT scalar(a[idx*NbPrimesQ+p]);

                // 
                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_0);
                tc->scalar.emplace_back(scalar * FieldT(b[7*(idx*NbPrimesQ)+0*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_1);
                tc->scalar.emplace_back(scalar * FieldT(b[7*(idx*NbPrimesQ)+1*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_2);
                tc->scalar.emplace_back(scalar * FieldT(b[7*(idx*NbPrimesQ)+2*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_3);
                tc->scalar.emplace_back(scalar * FieldT(b[7*(idx*NbPrimesQ)+3*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_4);
                tc->scalar.emplace_back(scalar * FieldT(b[7*(idx*NbPrimesQ)+4*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_5);
                tc->scalar.emplace_back(scalar * FieldT(b[7*(idx*NbPrimesQ)+5*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(blocOneIdx);
                tc->source_position.emplace_back(0);
                tc->scalar.emplace_back(scalar * FieldT(_publicData.special ? FieldT(3) : FieldT(0))
                                * FieldT(b[7*(idx*NbPrimesQ)+6*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_y_idx + idx);
                tc->scalar.emplace_back(scalar * minusOne * FieldT(c[(idx*NbPrimesQ) +p]));
                tc->target_position.emplace_back(blocCol);
            }

            {
                FieldT scalar(a_2[idx*NbPrimesQ +p]);

                // 
                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_2_0);
                tc->scalar.emplace_back(scalar * FieldT(b[7*(idx*NbPrimesQ)+0*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_2_1);
                tc->scalar.emplace_back(scalar * FieldT(b[7*(idx*NbPrimesQ)+1*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_2_2);
                tc->scalar.emplace_back(scalar * FieldT(b[7*(idx*NbPrimesQ)+2*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_2_3);
                tc->scalar.emplace_back(scalar * FieldT(b[7*(idx*NbPrimesQ)+3*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_2_4);
                tc->scalar.emplace_back(scalar * FieldT(b[7*(idx*NbPrimesQ)+4*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_2_5);
                tc->scalar.emplace_back(scalar * FieldT(b[7*(idx*NbPrimesQ)+5*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(blocOneIdx);
                tc->source_position.emplace_back(0);
                tc->scalar.emplace_back(scalar * FieldT(_publicData.special ? FieldT(3) : FieldT(0))
                                * FieldT(b[7*(idx*NbPrimesQ)+6*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_y_2_idx + idx);
                tc->scalar.emplace_back(scalar * minusOne * FieldT(c_2[(idx*NbPrimesQ) +p]));
                tc->target_position.emplace_back(blocCol);
            }

            {
                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_3_idx + idx);
                tc->scalar.emplace_back(FieldT(1));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(blocOneIdx);
                tc->source_position.emplace_back(0);
                tc->scalar.emplace_back(minusOne * a_3[(idx*NbPrimesQ)+p]);
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_y_3_idx + idx);
                tc->scalar.emplace_back(minusOne * FieldT(b_3[(idx*NbPrimesQ)+p]));
                tc->target_position.emplace_back(blocCol);
            }

            blocCol++;
        }
    }

    void addEquation11(
        transformation_constraint<FieldT> *tc, 
        size_t blocZeroIdx,
        size_t blocOneIdx,

        std::vector<uint64_t>& x_0,
        std::vector<uint64_t>& a_0,

        std::vector<uint64_t>& x_1,
        std::vector<uint64_t>& y_1,
        std::vector<uint64_t>& a_1,
        std::vector<uint64_t>& b_1,

        std::vector<uint64_t>& x_2,
        std::vector<uint64_t>& y_2,
        std::vector<uint64_t>& z_2,
        std::vector<uint64_t>& a_2,
        std::vector<uint64_t>& b_2,
        std::vector<uint64_t>& c_2,

        std::vector<uint64_t>& x_3,
        std::vector<uint64_t>& y_3,
        std::vector<uint64_t>& z_3,
        std::vector<uint64_t>& a_3,
        std::vector<uint64_t>& b_3,

        size_t& blocCol
    ) {
        // populating the misc bloc

        size_t& p = _publicData.modulusIdx;

        // per prime
        std::vector<FieldT> miscBloc(this->_blockSize, FieldT(0));
        size_t miscIdx=0;

        auto addEmptyVar = [&]() -> size_t {
            miscBloc[miscIdx++] = FieldT(0);
            return (miscIdx-1);
        };

        auto addEmptyVect = [&](std::vector<size_t> idxs) -> size_t {
            for (size_t idx : idxs) {miscBloc[miscIdx++] = FieldT(0);}
            return (miscIdx - idxs.size());
        };

        auto addVar = [&](uint64_t a) -> size_t {
            miscBloc[miscIdx++] = a;
            return (miscIdx-1);
        };

        auto addVect = [&](std::vector<uint64_t>& a, std::vector<size_t> idxs) -> size_t {
            for (size_t idx : idxs) {miscBloc[miscIdx++] = a[idx];}
            return (miscIdx - idxs.size());
        };

        constexpr int NbPrimesQ             = 21;
        constexpr int alphaCanLength        = 6;
        constexpr int alphaPSLength         = 6;
        constexpr int alphaGCDLength        = 18;

        // Declaring all variables
        size_t s_x_1_0;
        size_t s_x_1_1;
        size_t s_x_1_2;
        size_t s_x_1_3;
        size_t s_x_1_4;
        size_t s_x_1_5;

        size_t s_x_2_0;
        size_t s_x_2_1;
        size_t s_x_2_2;
        size_t s_x_2_3;
        size_t s_x_2_4;
        size_t s_x_2_5;

        size_t s_x_0_idx;

        size_t s_y_1_idx;

        size_t s_y_2_idx;
        size_t s_z_2_idx;

        size_t s_x_3_idx;
        size_t s_y_3_idx;
        size_t s_z_3_idx;

        // Vectors
        std::vector<size_t> v;
        for (size_t idx = 0; idx < alphaGCDLength; idx++) {v.push_back(idx*NbPrimesQ+p);}

        if (_secretData.isAvailable) {

            // Single Variables
            s_x_1_0   = addVar(x_1[              p]);
            s_x_1_1   = addVar(x_1[  NbPrimesQ + p]);
            s_x_1_2   = addVar(x_1[2*NbPrimesQ + p]);
            s_x_1_3   = addVar(x_1[3*NbPrimesQ + p]);
            s_x_1_4   = addVar(x_1[4*NbPrimesQ + p]);
            s_x_1_5   = addVar(x_1[5*NbPrimesQ + p]);

            s_x_2_0   = addVar(x_2[              p]);
            s_x_2_1   = addVar(x_2[  NbPrimesQ + p]);
            s_x_2_2   = addVar(x_2[2*NbPrimesQ + p]);
            s_x_2_3   = addVar(x_2[3*NbPrimesQ + p]);
            s_x_2_4   = addVar(x_2[4*NbPrimesQ + p]);
            s_x_2_5   = addVar(x_2[5*NbPrimesQ + p]);

            s_x_0_idx    = addVect(x_0, v);
            s_x_3_idx    = addVect(x_3, v);

            s_y_1_idx    = addVect(y_1, v);
            s_y_2_idx    = addVect(y_2, v);
            s_y_3_idx    = addVect(y_3, v);

            s_z_2_idx    = addVect(z_2, v);
            s_z_3_idx    = addVect(z_3, v);

        } else {
            // Single Variables
            s_x_1_0   = addEmptyVar();
            s_x_1_1   = addEmptyVar();
            s_x_1_2   = addEmptyVar();
            s_x_1_3   = addEmptyVar();
            s_x_1_4   = addEmptyVar();
            s_x_1_5   = addEmptyVar();

            s_x_2_0   = addEmptyVar();
            s_x_2_1   = addEmptyVar();
            s_x_2_2   = addEmptyVar();
            s_x_2_3   = addEmptyVar();
            s_x_2_4   = addEmptyVar();
            s_x_2_5   = addEmptyVar();

            s_x_0_idx    = addEmptyVect(v);
            s_x_3_idx    = addEmptyVect(v);

            s_y_1_idx    = addEmptyVect(v);
            s_y_2_idx    = addEmptyVect(v);
            s_y_3_idx    = addEmptyVect(v);

            s_z_2_idx    = addEmptyVect(v);
            s_z_3_idx    = addEmptyVect(v);
        }

        // Book-keep location of x shares
        _xsharesPSidxs_3.cols.assign({s_x_1_0, s_x_1_1, s_x_1_2, s_x_1_3, s_x_1_4, s_x_1_5});
        _ysharesPSidxs_3.cols.assign({s_x_2_0, s_x_2_1, s_x_2_2, s_x_2_3, s_x_2_4, s_x_2_5});

        // Populating inputs 
        size_t miscBlocIdx =_builder.add(new variable_constraint<FieldT>());
        _variables.add(miscBlocIdx,miscBloc);

        // Book-keeping for connecting proofs
        _xsharesPSidxs_3.blocs.assign(alphaPSLength, miscBlocIdx);
        _ysharesPSidxs_3.blocs.assign(alphaPSLength, miscBlocIdx);

        // xsharesGCD
        for (size_t idx = 0; idx < alphaGCDLength; idx++) this->_xsharesGCDidxs.cols.emplace_back(s_x_3_idx + idx);
        this->_xsharesGCDidxs.blocs.assign(alphaGCDLength, miscBlocIdx);

        // ysharesGCD
        for (size_t idx = 0; idx < alphaGCDLength; idx++) this->_ysharesGCDidxs.cols.emplace_back(s_x_0_idx + idx);
        this->_ysharesGCDidxs.blocs.assign(alphaGCDLength, miscBlocIdx);

        FieldT minusOne(FieldT::getModulus()-1);
        // Linear Combination
        for (size_t idx = 0; idx < alphaCanLength; idx++) {

            size_t i = idx * NbPrimesQ;

            // 11 - This constraint is broken down in three components for lisibility
            {
                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_0_idx + idx);
                tc->scalar.emplace_back(1);
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(blocOneIdx);
                tc->source_position.emplace_back(0);
                tc->scalar.emplace_back(a_0[i+p]);
                tc->target_position.emplace_back(blocCol);
            }

            {
                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_1_0);
                tc->scalar.emplace_back(minusOne * FieldT(a_1[7*i+0*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_1_1);
                tc->scalar.emplace_back(minusOne * FieldT(a_1[7*i+1*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_1_2);
                tc->scalar.emplace_back(minusOne * FieldT(a_1[7*i+2*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_1_3);
                tc->scalar.emplace_back(minusOne * FieldT(a_1[7*i+3*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_1_4);
                tc->scalar.emplace_back(minusOne * FieldT(a_1[7*i+4*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_1_5);
                tc->scalar.emplace_back(minusOne * FieldT(a_1[7*i+5*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(blocOneIdx);
                tc->source_position.emplace_back(0);
                tc->scalar.emplace_back(minusOne * FieldT(_publicData.special ? FieldT(3) : FieldT(0)) * FieldT(a_1[7*i+6*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_y_1_idx + idx);
                tc->scalar.emplace_back(b_1[i+p]);
                tc->target_position.emplace_back(blocCol);
            }

            {
                // 
                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_2_0);
                tc->scalar.emplace_back(minusOne * FieldT(a_2[7*i+0*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_2_1);
                tc->scalar.emplace_back(minusOne * FieldT(a_2[7*i+1*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_2_2);
                tc->scalar.emplace_back(minusOne * FieldT(a_2[7*i+2*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_2_3);
                tc->scalar.emplace_back(minusOne * FieldT(a_2[7*i+3*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_2_4);
                tc->scalar.emplace_back(minusOne * FieldT(a_2[7*i+4*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_2_5);
                tc->scalar.emplace_back(minusOne * FieldT(a_2[7*i+5*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(blocOneIdx);
                tc->source_position.emplace_back(0);
                tc->scalar.emplace_back(minusOne * FieldT(_publicData.special ? FieldT(3) : FieldT(0)) * FieldT(a_2[7*i+6*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_y_2_idx + idx);
                tc->scalar.emplace_back(b_2[i+p]);
                tc->target_position.emplace_back(blocCol);
            }

            {
                tc->block.emplace_back(blocOneIdx);
                tc->source_position.emplace_back(0);
                tc->scalar.emplace_back(_publicData.special ? 1 : 0);
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_z_2_idx + idx);
                tc->scalar.emplace_back(minusOne * FieldT(c_2[i+p]));
                tc->target_position.emplace_back(blocCol);
            }

            blocCol++;

            // 12- This is a different constraint, we need to ensure that both get evaluated to 0
            {
                tc->block.emplace_back(blocOneIdx);
                tc->source_position.emplace_back(0);
                tc->scalar.emplace_back(FieldT(a_3[i+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_3_idx + idx);
                tc->scalar.emplace_back(1);
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_y_3_idx + idx);
                tc->scalar.emplace_back(minusOne);
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_z_3_idx + idx);
                tc->scalar.emplace_back(minusOne * FieldT(b_3[i+p]));
                tc->target_position.emplace_back(blocCol);
            }

            blocCol++;
        }
    }

    void addEquation12(
        transformation_constraint<FieldT> *tc, 
        size_t blocZeroIdx,
        size_t blocOneIdx,

        std::vector<uint64_t>& a_0,

        std::vector<uint64_t>& x_1,
        std::vector<uint64_t>& y_1,
        std::vector<uint64_t>& z_1,
        std::vector<uint64_t>& t_1,
        std::vector<uint64_t>& a_1,
        std::vector<uint64_t>& b_1,
        std::vector<uint64_t>& c_1,

        std::vector<uint64_t>& x_2,
        std::vector<uint64_t>& y_2,
        std::vector<uint64_t>& z_2,
        std::vector<uint64_t>& a_2,
        std::vector<uint64_t>& b_2,

        size_t& blocCol
    ) {
        // populating the misc bloc

        size_t& p = _publicData.modulusIdx;

        // per prime
        std::vector<FieldT> miscBloc(this->_blockSize, FieldT(0));
        size_t miscIdx=0;

        auto addEmptyVar = [&]() -> size_t {
            miscBloc[miscIdx++] = FieldT(0);
            return (miscIdx-1);
        };

        auto addEmptyVect = [&](std::vector<size_t> idxs) -> size_t {
            for (size_t idx : idxs) {miscBloc[miscIdx++] = FieldT(0);}
            return (miscIdx - idxs.size());
        };

        auto addVar = [&](uint64_t a) -> size_t {
            miscBloc[miscIdx++] = a;
            return (miscIdx-1);
        };

        auto addVect = [&](std::vector<uint64_t>& a, std::vector<size_t> idxs) -> size_t {
            for (size_t idx : idxs) {miscBloc[miscIdx++] = a[idx];}
            return (miscIdx - idxs.size());
        };

        constexpr int NbPrimesQ             = 21;
        constexpr int alphaCanLength        = 6;
        constexpr int alphaPSLength         = 6;
        constexpr int alphaGCDLength        = 18;

        // Declaring all variables
        size_t s_x_1_0;
        size_t s_x_1_1;
        size_t s_x_1_2;
        size_t s_x_1_3;
        size_t s_x_1_4;
        size_t s_x_1_5;

        size_t s_z_1_0;
        size_t s_z_1_1;
        size_t s_z_1_2;
        size_t s_z_1_3;
        size_t s_z_1_4;
        size_t s_z_1_5;

        size_t s_y_1_idx;
        size_t s_t_1_idx;

        size_t s_x_2_idx;
        size_t s_y_2_idx;
        size_t s_z_2_idx;

        // Vectors
        std::vector<size_t> v;
        for (size_t idx = 0; idx < alphaGCDLength; idx++) {v.push_back(idx*NbPrimesQ+p);}

        if (_secretData.isAvailable) {

            // Single Variables
            s_x_1_0   = addVar(x_1[              p]);
            s_x_1_1   = addVar(x_1[  NbPrimesQ + p]);
            s_x_1_2   = addVar(x_1[2*NbPrimesQ + p]);
            s_x_1_3   = addVar(x_1[3*NbPrimesQ + p]);
            s_x_1_4   = addVar(x_1[4*NbPrimesQ + p]);
            s_x_1_5   = addVar(x_1[5*NbPrimesQ + p]);

            s_z_1_0   = addVar(z_1[              p]);
            s_z_1_1   = addVar(z_1[  NbPrimesQ + p]);
            s_z_1_2   = addVar(z_1[2*NbPrimesQ + p]);
            s_z_1_3   = addVar(z_1[3*NbPrimesQ + p]);
            s_z_1_4   = addVar(z_1[4*NbPrimesQ + p]);
            s_z_1_5   = addVar(z_1[5*NbPrimesQ + p]);

            s_y_1_idx    = addVect(y_1, v);
            s_t_1_idx    = addVect(t_1, v);

            s_x_2_idx    = addVect(x_2, v);
            s_y_2_idx    = addVect(y_2, v);
            s_z_2_idx    = addVect(z_2, v);

        } else {
            // Single Variables
            s_x_1_0   = addEmptyVar();
            s_x_1_1   = addEmptyVar();
            s_x_1_2   = addEmptyVar();
            s_x_1_3   = addEmptyVar();
            s_x_1_4   = addEmptyVar();
            s_x_1_5   = addEmptyVar();

            s_z_1_0   = addEmptyVar();
            s_z_1_1   = addEmptyVar();
            s_z_1_2   = addEmptyVar();
            s_z_1_3   = addEmptyVar();
            s_z_1_4   = addEmptyVar();
            s_z_1_5   = addEmptyVar();

            s_y_1_idx    = addEmptyVect(v);
            s_t_1_idx    = addEmptyVect(v);

            s_x_2_idx    = addEmptyVect(v);
            s_y_2_idx    = addEmptyVect(v);
            s_z_2_idx    = addEmptyVect(v);
        }

        // Book-keep location of x shares
        _xsharesPSidxs_4.cols.assign({s_x_1_0, s_x_1_1, s_x_1_2, s_x_1_3, s_x_1_4, s_x_1_5});
        _ysharesPSidxs_4.cols.assign({s_z_1_0, s_z_1_1, s_z_1_2, s_z_1_3, s_z_1_4, s_z_1_5});

        // Populating inputs 
        size_t miscBlocIdx =_builder.add(new variable_constraint<FieldT>());
        _variables.add(miscBlocIdx,miscBloc);

        // Book-keeping for connecting proofs
        _xsharesPSidxs_4.blocs.assign(alphaPSLength, miscBlocIdx);
        _ysharesPSidxs_4.blocs.assign(alphaPSLength, miscBlocIdx);

        FieldT minusOne(FieldT::getModulus()-1);

        // Linear Combination
        for (size_t idx = 0; idx < alphaCanLength; idx++) {

            size_t i = idx * NbPrimesQ;

            {
                tc->block.emplace_back(blocOneIdx);
                tc->source_position.emplace_back(0);
                tc->scalar.emplace_back(a_0[i+p]);
                tc->target_position.emplace_back(blocCol);
            }

            {
                FieldT d(a_1[i+p]);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_1_0);
                tc->scalar.emplace_back(minusOne * d * FieldT(b_1[7*i+0*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_1_1);
                tc->scalar.emplace_back(minusOne * d * FieldT(b_1[7*i+1*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_1_2);
                tc->scalar.emplace_back(minusOne * d * FieldT(b_1[7*i+2*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_1_3);
                tc->scalar.emplace_back(minusOne * d * FieldT(b_1[7*i+3*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_1_4);
                tc->scalar.emplace_back(minusOne * d * FieldT(b_1[7*i+4*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_1_5);
                tc->scalar.emplace_back(minusOne * d * FieldT(b_1[7*i+5*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(blocOneIdx);
                tc->source_position.emplace_back(0);
                tc->scalar.emplace_back(minusOne * d * FieldT(_publicData.special ? FieldT(3) : FieldT(0)) * FieldT(b_1[7*i+6*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_y_1_idx + idx);
                tc->scalar.emplace_back(d * FieldT(c_1[i + p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_z_1_0);
                tc->scalar.emplace_back(minusOne * d * FieldT(b_1[7*i+0*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_z_1_1);
                tc->scalar.emplace_back(minusOne * d * FieldT(b_1[7*i+1*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_z_1_2);
                tc->scalar.emplace_back(minusOne * d * FieldT(b_1[7*i+2*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_z_1_3);
                tc->scalar.emplace_back(minusOne * d * FieldT(b_1[7*i+3*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_z_1_4);
                tc->scalar.emplace_back(minusOne * d * FieldT(b_1[7*i+4*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_z_1_5);
                tc->scalar.emplace_back(minusOne * d * FieldT(b_1[7*i+5*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(blocOneIdx);
                tc->source_position.emplace_back(0);
                tc->scalar.emplace_back(minusOne * d * FieldT(_publicData.special ? FieldT(3) : FieldT(0)) * FieldT(b_1[7*i+6*NbPrimesQ+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_t_1_idx + idx);
                tc->scalar.emplace_back(d * FieldT(c_1[i + p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(blocOneIdx);
                tc->source_position.emplace_back(0);
                tc->scalar.emplace_back(d * FieldT(_publicData.special ? FieldT(1) : FieldT(0)));
                tc->target_position.emplace_back(blocCol);
            }

            {
                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_x_2_idx + idx);
                tc->scalar.emplace_back(minusOne * FieldT(a_2[i+p]));
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_y_2_idx + idx);
                tc->scalar.emplace_back(minusOne);
                tc->target_position.emplace_back(blocCol);

                tc->block.emplace_back(miscBlocIdx);
                tc->source_position.emplace_back(s_z_2_idx + idx);
                tc->scalar.emplace_back(FieldT(b_2[i+p]));
                tc->target_position.emplace_back(blocCol);
            }

            blocCol++;
        }
    }

    // Here we are adding constraints for an equation of the type:
    // a = x * b + y in the ring Zn/(x^n + 1)
    // or: q exists s.t. x * b + y - a - (x^n + 1) * q = 0

    // where    a is public
    //          b is public
    //          x is private
    //          y is private
    //          q is private

    // Returns the bloc indices for y, as these are often needed to tie in a separate constraint.

    void addEquation1(
            transformation_constraint<FieldT> *newTransformationConstraint, 
            std::vector<FieldT>& evaluationPoints, 
            size_t blocZeroIdx, 
            size_t blockOneIdx,
            std::vector<uint64_t> aIn,
            std::vector<uint64_t> bIn,
            std::vector<uint64_t> xIn,
            std::vector<uint64_t> yIn,
            std::vector<uint64_t> qIn,
            size_t& blocCol
            ) {

        // We assume that the bloc size divides the poly evaluation size for simplicity
        size_t nbBlocs = (_publicData.degree/_blockSize);    

        // Source All Data
        // ======================================================================================================

        // Public Data
        // A and bi are stored in the evaluation domain
        std::vector<FieldT> a = embed(aIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        std::vector<FieldT> b = embed(bIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);

        // Secret Data: either source or initialize if not available
        std::vector<FieldT> x;
        std::vector<FieldT> y;
        std::vector<FieldT> q;

        if (_secretData.isAvailable) {
            x = embed(xIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            y = embed(yIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            q = embed(qIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        } else {
            x.assign(_publicData.degree, FieldT(0));
            y.assign(_publicData.degree, FieldT(0));
            q.assign(_publicData.degree, FieldT(0));
        }

        // Building the constraint system
        // ======================================================================================================

        // Bounding Variables

        // Bounding Constraint
        // Rephrase the problem using the nearest 2^n
        // Here we look to bound by 10*sigma, so we need D larger than 2* (10*sigma) 
        // to phrase the problem as absolute numbers
        auto [d, C, D, log2D] = rephraseAs2Nct(_publicData.sigma*10);

        // Bounding variables
        vector<size_t> xCoefsBlocIds;
        if (_publicData.modulusIdx == 0) {
            xCoefsBlocIds = addBoundingConstraint(_publicData.degree, x, FieldT(C), FieldT(d), static_cast<uint64_t>(log2D), blockOneIdx);
        } else {
            xCoefsBlocIds = addSplitVariable(x);
        }

        // Splitting variables
        // vector<size_t> xCoefsBlocIds = addSplitVariable(x);
        vector<size_t> yCoefsBlocIds = addSplitVariable(y);
        vector<size_t> qCoefsBlocIds = addSplitVariable(q);

        // Evaluate on specific points
        DBG("Number of evaluation points:" << evaluationPoints.size());
        for (FieldT eval : evaluationPoints) {
            std::vector<FieldT> powersOfEval(_publicData.degree);

            for (size_t powIdx = 0; powIdx < _publicData.degree; powIdx++) {
                powersOfEval[powIdx] = eval^powIdx;
            }

            // Adding x * b + y - a = (x^n + 1) * q constraint 
            assert (blocCol < this->_blockSize);

            addEvalPt(newTransformationConstraint, blocCol, naiveEval(b, powersOfEval), xCoefsBlocIds, powersOfEval); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(1), yCoefsBlocIds, powersOfEval); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1)*naiveEval(a, powersOfEval), blockOneIdx); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1)*(eval*powersOfEval.back()+1), qCoefsBlocIds, powersOfEval); 

            blocCol++;
        }

    }

    // Here we are adding constraints for an equation of the type:
    // a = x * b + y + z in the ring Zn/(x^n + 1)
    // or: q exists s.t. x * b + y + z - a - (x^n + 1) * q = 0

    // where    a is public
    //          b is public
    //          x is private
    //          y is private
    //          z is private
    //          q is private

    // Returns the bloc indices for y, as these are often needed to tie in a separate constraint.

    void addEquation2(
            transformation_constraint<FieldT> *newTransformationConstraint, 
            std::vector<FieldT>& evaluationPoints, 
            size_t blocZeroIdx, 
            size_t blockOneIdx,
            std::vector<uint64_t> aIn,
            std::vector<uint64_t> bIn,
            std::vector<uint64_t> xIn,
            std::vector<uint64_t> yIn,
            std::vector<uint64_t> zIn, //xi
            std::vector<uint64_t> qIn,
            size_t& blocCol
            ) {

        // We assume that the bloc size divides the poly evaluation size for simplicity
        size_t nbBlocs = (_publicData.degree/_blockSize);    

        // Source All Data
        // ======================================================================================================

        // Public Data
        // A and bi are stored in the evaluation domain
        std::vector<FieldT> a = embed(aIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        std::vector<FieldT> b = embed(bIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);

        // Secret Data: either source or initialize if not available
        std::vector<FieldT> x;
        std::vector<FieldT> y;
        std::vector<FieldT> z;
        std::vector<FieldT> q;

        if (_secretData.isAvailable) {
            x = embed(xIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            y = embed(yIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            z = embed(zIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            q = embed(qIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        } else {
            x.assign(_publicData.degree, FieldT(0));
            y.assign(_publicData.degree, FieldT(0));
            z.assign(_publicData.degree, FieldT(0));
            q.assign(_publicData.degree, FieldT(0));
        }

        // Building the constraint system
        // ======================================================================================================

        // Splitting variables
        vector<size_t> xCoefsBlocIds = addSplitVariable(x);
        vector<size_t> yCoefsBlocIds = addSplitVariable(y);
        vector<size_t> zCoefsBlocIds = addSplitVariable(z);

        this->_xprimeIds.assign(zCoefsBlocIds.begin(), zCoefsBlocIds.end());
        
        vector<size_t> qCoefsBlocIds = addSplitVariable(q);

        // Evaluate on specific points
        DBG("Number of evaluation points:" << evaluationPoints.size());
        for (FieldT eval : evaluationPoints) {
            std::vector<FieldT> powersOfEval(_publicData.degree);

            for (size_t powIdx = 0; powIdx < _publicData.degree; powIdx++) {
                powersOfEval[powIdx] = eval^powIdx;
            }

            assert (blocCol < this->_blockSize);

            // Adding x * b + y + z - a = (x^n + 1) * q constraint 
            addEvalPt(newTransformationConstraint, blocCol, naiveEval(b, powersOfEval), xCoefsBlocIds, powersOfEval); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(1), yCoefsBlocIds, powersOfEval); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(1), zCoefsBlocIds, powersOfEval); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1)*naiveEval(a, powersOfEval), blockOneIdx); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1)*(eval*powersOfEval.back()+1), qCoefsBlocIds, powersOfEval); 

            blocCol++;
        }
    }

    // Here we are adding constraints for an equation of the type:
    // t = x * b + y + z in the ring Zn/(x^n + 1)
    // or: q exists s.t. x * b + y + z - t - (x^n + 1) * q = 0

    // where    a is public
    //          b is public
    //          x is private
    //          y is private
    //          z is private
    //          q is private

    std::vector<size_t> addEquation2p(
            transformation_constraint<FieldT> *newTransformationConstraint, 
            std::vector<FieldT>& evaluationPoints, 
            size_t blocZeroIdx, 
            size_t blockOneIdx,
            std::vector<size_t> tCoefsBlocIds,
            std::vector<uint64_t> bIn,
            std::vector<uint64_t> xIn,
            std::vector<uint64_t> yIn,
            std::vector<uint64_t> zIn,
            std::vector<uint64_t> qIn,
            size_t& blocCol
            ) {

        // We assume that the bloc size divides the poly evaluation size for simplicity
        size_t nbBlocs = (_publicData.degree/_blockSize);    

        // Source All Data
        // ======================================================================================================

        // Public Data
        // A and bi are stored in the evaluation domain
        std::vector<FieldT> b = embed(bIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);

        // Secret Data: either source or initialize if not available
        std::vector<FieldT> x;
        std::vector<FieldT> y;
        std::vector<FieldT> z;
        std::vector<FieldT> q;

        if (_secretData.isAvailable) {
            x = embed(xIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            y = embed(yIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            z = embed(zIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            q = embed(qIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        } else {
            x.assign(_publicData.degree, FieldT(0));
            y.assign(_publicData.degree, FieldT(0));
            z.assign(_publicData.degree, FieldT(0));
            q.assign(_publicData.degree, FieldT(0));
        }

        // Building the constraint system
        // ======================================================================================================

        // Splitting variables
        vector<size_t> xCoefsBlocIds = addSplitVariable(x);
        vector<size_t> yCoefsBlocIds = addSplitVariable(y);
        vector<size_t> zCoefsBlocIds = addSplitVariable(z);
        vector<size_t> qCoefsBlocIds = addSplitVariable(q);

        // Evaluate on specific points
        DBG("Number of evaluation points:" << evaluationPoints.size());
        for (FieldT eval : evaluationPoints) {
            std::vector<FieldT> powersOfEval(_publicData.degree);

            for (size_t powIdx = 0; powIdx < _publicData.degree; powIdx++) {
                powersOfEval[powIdx] = eval^powIdx;
            }

            assert (blocCol < this->_blockSize);

            // Adding x * b + y + z - a = (x^n + 1) * q constraint 
            addEvalPt(newTransformationConstraint, blocCol, naiveEval(b, powersOfEval), xCoefsBlocIds, powersOfEval); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(1), yCoefsBlocIds, powersOfEval); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(1), zCoefsBlocIds, powersOfEval); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1), tCoefsBlocIds, powersOfEval); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1)*(eval*powersOfEval.back()+1), qCoefsBlocIds, powersOfEval); 

            blocCol++;
        }
    }

    // Here we are adding constraints for an equation of the type:
    // a = x * b - y * c - z in the ring Zn/(x^n + 1)
    // or: q exists s.t. x * b - y * c - z - a - (x^n + 1) * q = 0

    // where    a is public
    //          b is public
    //          c is public

    //          x is private
    //          y is private
    //          z is private

    //          q is private

    // Returns the bloc indices for y, as these are often needed to tie in a separate constraint.

    void addEquation3(
            transformation_constraint<FieldT> *newTransformationConstraint, 
            std::vector<FieldT>& evaluationPoints, 
            size_t blocZeroIdx, 
            size_t blockOneIdx,
            std::vector<uint64_t> aIn,
            std::vector<uint64_t> bIn,
            std::vector<uint64_t> cIn,
            std::vector<uint64_t> xIn,
            std::vector<uint64_t> yIn,
            std::vector<uint64_t> zIn,
            std::vector<uint64_t> qIn,
            size_t& blocCol
            ) {

        // We assume that the bloc size divides the poly evaluation size for simplicity
        size_t nbBlocs = (_publicData.degree/_blockSize);    

        // Source All Data
        // ======================================================================================================

        // Public Data
        // A and bi are stored in the evaluation domain
        std::vector<FieldT> a = embed(aIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        std::vector<FieldT> b = embed(bIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        std::vector<FieldT> c = embed(cIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);

        // Secret Data: either source or initialize if not available
        std::vector<FieldT> x;
        std::vector<FieldT> y;
        std::vector<FieldT> z;
        std::vector<FieldT> q;

        if (_secretData.isAvailable) {
            x = embed(xIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            y = embed(yIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            z = embed(zIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            q = embed(qIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        } else {
            x.assign(_publicData.degree, FieldT(0));
            y.assign(_publicData.degree, FieldT(0));
            z.assign(_publicData.degree, FieldT(0));
            q.assign(_publicData.degree, FieldT(0));
        }

        // Building the constraint system
        // ======================================================================================================

        // Bounding Variables

        // Bounding Constraint
        // Rephrase the problem using the nearest 2^n
        // Here we look to bound by 10*sigma, so we need D larger than 2* (10*sigma) 
        // to phrase the problem as absolute numbers
        auto [d, C, D, log2D] = rephraseAs2Nct(_publicData.sigma*10);

        // Bounding variables
        vector<size_t> yCoefsBlocIds;
        if (_publicData.modulusIdx == 0) {
            yCoefsBlocIds = addBoundingConstraint(_publicData.degree, y, FieldT(C), FieldT(d), static_cast<uint64_t>(log2D), blockOneIdx);
        } else {
            yCoefsBlocIds = addSplitVariable(y);
        }

        // Splitting variables
        vector<size_t> xCoefsBlocIds = addSplitVariable(x);
        vector<size_t> zCoefsBlocIds = addSplitVariable(z);
        vector<size_t> qCoefsBlocIds = addSplitVariable(q);

        this->_yprimeIds.assign(xCoefsBlocIds.begin(), xCoefsBlocIds.end());

        // Evaluate on specific points
        DBG("Number of evaluation points:" << evaluationPoints.size());
        for (FieldT eval : evaluationPoints) {
            std::vector<FieldT> powersOfEval(_publicData.degree);

            for (size_t powIdx = 0; powIdx < _publicData.degree; powIdx++) {
                powersOfEval[powIdx] = eval^powIdx;
            }

            assert (blocCol < this->_blockSize);

            addEvalPt(newTransformationConstraint, blocCol, naiveEval(b, powersOfEval), xCoefsBlocIds, powersOfEval);
            addEvalPt(newTransformationConstraint, blocCol, FieldT(0) - naiveEval(c, powersOfEval), yCoefsBlocIds, powersOfEval);
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1), zCoefsBlocIds, powersOfEval);
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1)*naiveEval(a, powersOfEval), blockOneIdx);
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1)*(eval*powersOfEval.back()+1), qCoefsBlocIds, powersOfEval);

            blocCol++;
        }
    }

    // Here we are adding constraints for an equation of the type:
    // a = x * b + z - y * c - t - m in the ring Zn/(x^n + 1)
    // or: q exists s.t. x * b + z - y * c - t - m - (x^n + 1) * q = 0

    // where    a is public
    //          b is public
    //          c is public

    //          x is private
    //          y is private
    //          z is private
    //          t is private
    //          m is private

    //          q is private

    // Returns the bloc indices for y, as these are often needed to tie in a separate constraint.

    void addEquation4(
            transformation_constraint<FieldT> *newTransformationConstraint, 
            std::vector<FieldT>& evaluationPoints, 
            size_t blocZeroIdx, 
            size_t blockOneIdx,
            std::vector<uint64_t> aIn,
            std::vector<uint64_t> bIn,
            std::vector<uint64_t> cIn,
            std::vector<uint64_t> xIn,
            std::vector<uint64_t> yIn,
            std::vector<uint64_t> zIn,
            std::vector<uint64_t> tIn,
            std::vector<uint64_t> mIn,
            std::vector<uint64_t> qIn,
            size_t& blocCol
            ) {

        // We assume that the bloc size divides the poly evaluation size for simplicity
        size_t nbBlocs = (_publicData.degree/_blockSize);    

        // Source All Data
        // ======================================================================================================

        // Public Data
        // A and bi are stored in the evaluation domain
        std::vector<FieldT> a = embed(aIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        std::vector<FieldT> b = embed(bIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        std::vector<FieldT> c = embed(cIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);

        // Secret Data: either source or initialize if not available
        std::vector<FieldT> x;
        std::vector<FieldT> y;
        std::vector<FieldT> z;
        std::vector<FieldT> t;
        std::vector<FieldT> m;
        std::vector<FieldT> q;

        if (_secretData.isAvailable) {
            x = embed(xIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            y = embed(yIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            z = embed(zIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            t = embed(tIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            m = embed(mIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            q = embed(qIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        } else {
            x.assign(_publicData.degree, FieldT(0));
            y.assign(_publicData.degree, FieldT(0));
            z.assign(_publicData.degree, FieldT(0));
            t.assign(_publicData.degree, FieldT(0));
            m.assign(_publicData.degree, FieldT(0));
            q.assign(_publicData.degree, FieldT(0));
        }

        // Building the constraint system
        // ======================================================================================================

        // Splitting variables
        vector<size_t> xCoefsBlocIds = addSplitVariable(x);
        vector<size_t> yCoefsBlocIds = addSplitVariable(y);
        vector<size_t> zCoefsBlocIds = addSplitVariable(z);
        vector<size_t> tCoefsBlocIds = addSplitVariable(t);
        vector<size_t> mCoefsBlocIds = addSplitVariable(m);
        vector<size_t> qCoefsBlocIds = addSplitVariable(q);

        this->_zprimeIds.assign(mCoefsBlocIds.begin(), mCoefsBlocIds.end());

        // Evaluate on specific points
        DBG("Number of evaluation points:" << evaluationPoints.size());
        for (FieldT eval : evaluationPoints) {
            std::vector<FieldT> powersOfEval(_publicData.degree);

            for (size_t powIdx = 0; powIdx < _publicData.degree; powIdx++) {
                powersOfEval[powIdx] = eval^powIdx;
            }

            assert (blocCol < this->_blockSize);

            addEvalPt(newTransformationConstraint, blocCol, naiveEval(b, powersOfEval), xCoefsBlocIds, powersOfEval);
            addEvalPt(newTransformationConstraint, blocCol, FieldT(1), zCoefsBlocIds, powersOfEval);
            addEvalPt(newTransformationConstraint, blocCol, FieldT(0) - naiveEval(c, powersOfEval), yCoefsBlocIds, powersOfEval);
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1), tCoefsBlocIds, powersOfEval);
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1), mCoefsBlocIds, powersOfEval);
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1)*naiveEval(a, powersOfEval), blockOneIdx);
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1)*(eval*powersOfEval.back()+1), qCoefsBlocIds, powersOfEval);

            blocCol++;
        }
    }

    // Here we are adding constraints for an equation of the type:
    // a = c - x * b  + y in the ring Zn/(x^n + 1)
    // or: q exists s.t. c - x * b  + y - a - (x^n + 1) * q = 0

    // where    a is public
    //          b is public
    //          c is public
    //          x is private
    //          y is private
    //          q is private

    // Returns the bloc indices for y, as these are often needed to tie in a separate constraint.

    void addEquation5(
            transformation_constraint<FieldT> *newTransformationConstraint, 
            std::vector<FieldT>& evaluationPoints, 
            size_t blocZeroIdx, 
            size_t blockOneIdx,
            std::vector<uint64_t> aIn,
            std::vector<uint64_t> bIn,
            std::vector<uint64_t> cIn,
            std::vector<uint64_t> xIn,
            std::vector<uint64_t> yIn,
            std::vector<uint64_t> qIn,
            size_t& blocCol
            ) {

        // We assume that the bloc size divides the poly evaluation size for simplicity
        size_t nbBlocs = (_publicData.degree/_blockSize);    

        // Source All Data
        // ======================================================================================================

        // Public Data
        // A and bi are stored in the evaluation domain
        std::vector<FieldT> a = embed(aIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        std::vector<FieldT> b = embed(bIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        std::vector<FieldT> c = embed(cIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);

        // Secret Data: either source or initialize if not available
        std::vector<FieldT> x;
        std::vector<FieldT> y;
        std::vector<FieldT> q;

        if (_secretData.isAvailable) {
            x = embed(xIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            y = embed(yIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
            q = embed(qIn, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
        } else {
            x.assign(_publicData.degree, FieldT(0));
            y.assign(_publicData.degree, FieldT(0));
            q.assign(_publicData.degree, FieldT(0));
        }

        // Building the constraint system
        // ======================================================================================================

        // Splitting variables
        vector<size_t> xCoefsBlocIds = addSplitVariable(x);
        vector<size_t> yCoefsBlocIds = addSplitVariable(y);
        vector<size_t> qCoefsBlocIds = addSplitVariable(q);

        // Evaluate on specific points
        DBG("Number of evaluation points:" << evaluationPoints.size());
        for (FieldT eval : evaluationPoints) {
            std::vector<FieldT> powersOfEval(_publicData.degree,FieldT(0));

            for (size_t powIdx = 0; powIdx < _publicData.degree; powIdx++) {
                powersOfEval[powIdx] = eval^powIdx;
            }

            assert (blocCol < this->_blockSize);

            if (_publicData.special)
            {
                addEvalPt(newTransformationConstraint, blocCol, naiveEval(c, powersOfEval), blockOneIdx); 
            }

            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1)*naiveEval(b, powersOfEval), xCoefsBlocIds, powersOfEval); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(1), yCoefsBlocIds, powersOfEval); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1)*naiveEval(a, powersOfEval), blockOneIdx); 
            addEvalPt(newTransformationConstraint, blocCol, FieldT(FieldT::getModulus()-1)*(eval*powersOfEval.back()+1), qCoefsBlocIds, powersOfEval); 

            blocCol++;
        }
    }

    void addEquationJacobi(
            transformation_constraint<FieldT> *newTransformationConstraint, 
            std::vector<FieldT>& evaluationPoints, 
            size_t blocZeroIdx, 
            size_t blockOneIdx,
            const mpz_class &z,
            const mpz_class &r,
            const mpz_class &x,
            const mpz_class &e,
            size_t& blocCol
            ) {

        // We assume that the bloc size divides the poly evaluation size for simplicity
        size_t nbBlocs = (_publicData.degree/_blockSize);    

        // Building the constraint system
        // ======================================================================================================

        // Splitting variables

        // Evaluate on specific points
        DBG("Number of evaluation points:" << evaluationPoints.size());
    }

    void addIntrablocBoundingConstraintLarge(std::vector<FieldT> varsIn, std::vector<FieldT> XplusCs, std::vector<FieldT> XplusDs, std::vector<FieldT> Cs, std::vector<FieldT> ds, std::vector<uint64_t> log2Ds, std::vector<std::vector<FieldT>> bitDecompositionsIn, size_t blockOneIdx, size_t blockZeroIdx) {
        std::vector<FieldT> vars;
        std::vector<size_t> bitDecompositionReferences;
        std::vector<FieldT> bitDecompositions;
        size_t colVar = 0;
        size_t colBit = 0;
        FieldT bitDecompositionScalar[MAX_BITS_SIZE];
        size_t initialOffset = varsIn.size();

        FieldT powerOfTwo = 1;
        for (size_t i = 0; i < MAX_BITS_SIZE; i++) {

            bitDecompositionScalar[i] = powerOfTwo;
            powerOfTwo = powerOfTwo * 2;
        }

        vars.assign(varsIn.begin(), varsIn.end());

        for (size_t index = 0; index < XplusCs.size(); index++) {

            // Create bits representation and retain index corresponding to their position; update appropriately based on bit length
            // Populate bitDecompositions vector

            // Note down the position for this bit representation 
            bitDecompositionReferences.emplace_back(bitDecompositions.size());

            // For more compact representation, we will add both representations within the same variable,
            // interleaving these inside blocs
            vars.emplace_back(XplusCs[index]);
            vars.emplace_back(XplusDs[index]);

            for (size_t bitIndex = 0; bitIndex < log2Ds[index]; bitIndex++) {

                // interleaving bit decompositions as well 
                bitDecompositions.emplace_back(bitDecompositionsIn[index*2][bitIndex]);
                bitDecompositions.emplace_back(bitDecompositionsIn[index*2+1][bitIndex]);
            }
        }

        // Populate the witness accordingly
        std::vector<size_t> varsIndices                = addSplitVariable(vars);
        std::vector<size_t> bitDecompositionsIndices   = addSplitVariable(bitDecompositions);

        // Enter the quadratic constraints for the bit blocs
        for (auto bitRepresentationBlocIdx : bitDecompositionsIndices) {

            // Quadratic constraints for bits
            // We don't need to create a variable for that, this constraint should equal 0
            size_t bitSquaredRepresentationIdx = _builder.add(new quadratic_constraint<FieldT>(
                        {bitRepresentationBlocIdx},
                        {bitRepresentationBlocIdx})
                        ); 

            _builder.ensure_equality(bitRepresentationBlocIdx, bitSquaredRepresentationIdx);
        }

        // Enter the intrabloc constraints for the decompositions
        //          browse through each varC/varD: 
        //          retrieve decomposition 

        auto writeIntrablocConstraint = [&](transformation_constraint<FieldT> *ct, size_t varsIndex, size_t actualTest, size_t col)->void {

            // Adding the portion of the constraint tied to individual bits
            for (size_t bitRepresentationIndex = 0; bitRepresentationIndex < log2Ds[varsIndex]; bitRepresentationIndex++) {
                size_t offsetIndex = bitDecompositionReferences[varsIndex] + (2*bitRepresentationIndex) + actualTest;

                ct->block.emplace_back(bitDecompositionsIndices[static_cast<size_t>(static_cast<float>(offsetIndex) / static_cast<float>(_blockSize))]);
                ct->scalar.emplace_back(bitDecompositionScalar[bitRepresentationIndex]);
                ct->source_position.emplace_back(offsetIndex % _blockSize);
                ct->target_position.emplace_back(col);
            }

            // Minus the actual fully formed number
            ct->block.emplace_back(varsIndices[static_cast<size_t>(static_cast<float>(initialOffset + (varsIndex*2) + actualTest) / static_cast<float>(_blockSize))]);
            ct->scalar.emplace_back(FieldT(FieldT::getModulus() - 1));
            ct->source_position.emplace_back((initialOffset + (varsIndex*2) + actualTest) % _blockSize);
            ct->target_position.emplace_back(col);
        };

        auto writeAdditionConstraint = [&](transformation_constraint<FieldT> *ct, size_t varsIndex, size_t actualTest, size_t col)->void {

            // X
            ct->block.emplace_back(varsIndices[static_cast<size_t>(static_cast<float>(varsIndex) / static_cast<float>(_blockSize))]);
            ct->scalar.emplace_back(FieldT(1));
            ct->source_position.emplace_back(varsIndex % _blockSize);
            ct->target_position.emplace_back(col);

            // C or d
            ct->block.emplace_back(blockOneIdx);
            ct->scalar.emplace_back(actualTest == 0 ? Cs[varsIndex] : ds[varsIndex]);
            ct->source_position.emplace_back(0);
            ct->target_position.emplace_back(col);

            // X + C or X + d
            ct->block.emplace_back(varsIndices[static_cast<size_t>(static_cast<float>(initialOffset + (varsIndex*2) + actualTest) / static_cast<float>(_blockSize))]);
            ct->scalar.emplace_back(FieldT(FieldT::getModulus() - 1));
            ct->source_position.emplace_back((initialOffset + (varsIndex*2) + actualTest) % _blockSize);
            ct->target_position.emplace_back(col);
        };

        transformation_constraint<FieldT> *ct;

        size_t col = 0;
        // Bit Decomposition Constraints on X + C
        for (size_t varsIndex = 0; varsIndex < varsIn.size(); varsIndex++) {
            if (col == 0) {ct = new transformation_constraint<FieldT>;}

            writeIntrablocConstraint(ct, varsIndex, 0, col++);

            if (col > _blockSize) {
                size_t bitRecompositionIdx = _builder.add(ct);
                _builder.ensure_equality(bitRecompositionIdx, blockZeroIdx);
                col = 0;
            }
        }

        // Bit Decomposition Constraint on X + d
        for (size_t varsIndex = 0; varsIndex < varsIn.size(); varsIndex++) {
            if (col == 0) {ct = new transformation_constraint<FieldT>;}

            writeIntrablocConstraint(ct, varsIndex, 1, col++);

            if (col > _blockSize) {
                size_t bitRecompositionIdx = _builder.add(ct);
                _builder.ensure_equality(bitRecompositionIdx, blockZeroIdx);
                col = 0;
            }
        }

        // Constraints on X + C
        for (size_t varsIndex = 0; varsIndex < varsIn.size(); varsIndex++) {
            if (col == 0) {ct = new transformation_constraint<FieldT>;}

            writeAdditionConstraint(ct, varsIndex, 0, col++);

            if (col > _blockSize) {
                size_t miscConstraintsIdx = _builder.add(ct);
                _builder.ensure_equality(miscConstraintsIdx, blockZeroIdx);
                col = 0;
            }
        }

        // Constraints on X + d
        for (size_t varsIndex = 0; varsIndex < varsIn.size(); varsIndex++) {
            if (col == 0) {ct = new transformation_constraint<FieldT>;}

            writeAdditionConstraint(ct, varsIndex, 1, col++);

            if (col > _blockSize) {
                size_t miscConstraintsIdx = _builder.add(ct);
                _builder.ensure_equality(miscConstraintsIdx, blockZeroIdx);
                col = 0;
            }
        }

        if (col != 0) {
                size_t miscConstraintsIdx = _builder.add(ct);
                // padBloc(ct, col, blockZeroIdx);

                _builder.ensure_equality(miscConstraintsIdx, blockZeroIdx);
            }
    }

    void addIntrablocBoundingConstraint(std::vector<FieldT> varsIn, std::vector<FieldT> ds, size_t blockOneIdx, size_t blockZeroIdx) {
        std::vector<FieldT> vars;
        std::vector<FieldT> Cs;
        std::vector<size_t> bitDecompositionReferences;
        std::vector<FieldT> bitDecompositions;
        std::vector<uint64_t> log2Ds;
        size_t colVar = 0;
        size_t colBit = 0;
        FieldT bitDecompositionScalar[MAX_BITS_SIZE];
        size_t initialOffset = varsIn.size();

        FieldT powerOfTwo = 1; 
        for (size_t i = 0; i < MAX_BITS_SIZE; i++) {

            bitDecompositionScalar[i] = powerOfTwo;
            powerOfTwo = powerOfTwo * 2;
        }

        vars.assign(varsIn.begin(), varsIn.end());

        for (auto index = 0; index < varsIn.size(); index++) {

            // Rephrase constraint, store Cs, Ds, log2Ds
            auto [d, C, D, log2D] = rephraseAs2Nct(static_cast<uint64_t>(ds[index].getValue()));

            // Create bits representation and retain index corresponding to their position; update appropriately based on bit length
            // Populate bitDecompositions vector

            // Note down the position for this bit representation 
            bitDecompositionReferences.emplace_back(bitDecompositions.size());
            log2Ds.emplace_back(log2D);
            Cs.emplace_back(C);

            // For more compact representation, we will add both representations within the same variable,
            // interleaving these inside blocs
            FieldT XplusC = vars[index] + C;
            FieldT Xplusd = vars[index] + d;

            vars.emplace_back(XplusC);
            vars.emplace_back(Xplusd);

            for (size_t bitIndex = 0; bitIndex < log2D; bitIndex++) {
                uint64_t bitPlusC = static_cast<uint64_t>((static_cast<uint64_t>(XplusC.getValue()) & (1ull << bitIndex)) >> bitIndex);
                uint64_t bitPlusD = static_cast<uint64_t>((static_cast<uint64_t>(Xplusd.getValue()) & (1ull << bitIndex)) >> bitIndex);

                // interleaving bit decompositions as well
                bitDecompositions.emplace_back(bitPlusC);
                bitDecompositions.emplace_back(bitPlusD);
            }
        }

        // Populate the witness accordingly
        std::vector<size_t> varsIndices                = addSplitVariable(vars);
        std::vector<size_t> bitDecompositionsIndices   = addSplitVariable(bitDecompositions);

        // Enter the quadratic constraints for the bit blocs
        for (auto bitRepresentationBlocIdx : bitDecompositionsIndices) {

            // Quadratic constraints for bits
            // We don't need to create a variable for that, this constraint should equal 0
            size_t bitSquaredRepresentationIdx = _builder.add(new quadratic_constraint<FieldT>(
                        {bitRepresentationBlocIdx},
                        {bitRepresentationBlocIdx})
                        ); 

            _builder.ensure_equality(bitRepresentationBlocIdx, bitSquaredRepresentationIdx);
        }

        // Enter the intrabloc constraints for the decompositions
        //          browse through each varC/varD: 
        //          retrieve decomposition 

        auto writeIntrablocConstraint = [&](transformation_constraint<FieldT> *ct, size_t varsIndex, size_t actualTest, size_t col)->void {

            // Adding the portion of the constraint tied to individual bits
            for (size_t bitRepresentationIndex = 0; bitRepresentationIndex < log2Ds[varsIndex]; bitRepresentationIndex++) {
                size_t offsetIndex = bitDecompositionReferences[varsIndex] + (2*bitRepresentationIndex) + actualTest;

                ct->block.emplace_back(bitDecompositionsIndices[static_cast<size_t>(static_cast<float>(offsetIndex) / static_cast<float>(_blockSize))]);
                ct->scalar.emplace_back(bitDecompositionScalar[bitRepresentationIndex]);
                ct->source_position.emplace_back(offsetIndex % _blockSize);
                ct->target_position.emplace_back(col);
            }

            // Minus the actual fully formed number
            ct->block.emplace_back(varsIndices[static_cast<size_t>(static_cast<float>(initialOffset + (varsIndex*2) + actualTest) / static_cast<float>(_blockSize))]);
            ct->scalar.emplace_back(FieldT(FieldT::getModulus() - 1));
            ct->source_position.emplace_back((initialOffset + (varsIndex*2) + actualTest) % _blockSize);
            ct->target_position.emplace_back(col);
        };

        auto writeAdditionConstraint = [&](transformation_constraint<FieldT> *ct, size_t varsIndex, size_t actualTest, size_t col)->void {

            // X
            ct->block.emplace_back(varsIndices[static_cast<size_t>(static_cast<float>(varsIndex) / static_cast<float>(_blockSize))]);
            ct->scalar.emplace_back(FieldT(1));
            ct->source_position.emplace_back(varsIndex % _blockSize);
            ct->target_position.emplace_back(col);

            // C or d
            ct->block.emplace_back(blockOneIdx);
            ct->scalar.emplace_back(actualTest == 0 ? Cs[varsIndex] : ds[varsIndex]);
            ct->source_position.emplace_back(0);
            ct->target_position.emplace_back(col);

            // X + C or X + d
            ct->block.emplace_back(varsIndices[static_cast<size_t>(static_cast<float>(initialOffset + (varsIndex*2) + actualTest) / static_cast<float>(_blockSize))]);
            ct->scalar.emplace_back(FieldT(FieldT::getModulus() - 1));
            ct->source_position.emplace_back((initialOffset + (varsIndex*2) + actualTest) % _blockSize);
            ct->target_position.emplace_back(col);
        };

        transformation_constraint<FieldT> *ct;

        size_t col = 0;
        // Bit Decomposition Constraints on X + C
        for (size_t varsIndex = 0; varsIndex < varsIn.size(); varsIndex++) {
            if (col == 0) {ct = new transformation_constraint<FieldT>;}

            writeIntrablocConstraint(ct, varsIndex, 0, col++);

            if (col > _blockSize) {
                size_t bitRecompositionIdx = _builder.add(ct);
                _builder.ensure_equality(bitRecompositionIdx, blockZeroIdx);
                col = 0;
            }
        }

        // Bit Decomposition Constraint on X + d
        for (size_t varsIndex = 0; varsIndex < varsIn.size(); varsIndex++) {
            if (col == 0) {ct = new transformation_constraint<FieldT>;}

            writeIntrablocConstraint(ct, varsIndex, 1, col++);

            if (col > _blockSize) {
                size_t bitRecompositionIdx = _builder.add(ct);
                _builder.ensure_equality(bitRecompositionIdx, blockZeroIdx);
                col = 0;
            }
        }

        // Constraints on X + C
        for (size_t varsIndex = 0; varsIndex < varsIn.size(); varsIndex++) {
            if (col == 0) {ct = new transformation_constraint<FieldT>;}

            writeAdditionConstraint(ct, varsIndex, 0, col++);

            if (col > _blockSize) {
                size_t miscConstraintsIdx = _builder.add(ct);
                _builder.ensure_equality(miscConstraintsIdx, blockZeroIdx);
                col = 0;
            }
        }

        // Constraints on X + d
        for (size_t varsIndex = 0; varsIndex < varsIn.size(); varsIndex++) {
            if (col == 0) {ct = new transformation_constraint<FieldT>;}

            writeAdditionConstraint(ct, varsIndex, 1, col++);

            if (col > _blockSize) {
                size_t miscConstraintsIdx = _builder.add(ct);
                _builder.ensure_equality(miscConstraintsIdx, blockZeroIdx);
                col = 0;
            }
        }

        if (col != 0) {
                size_t miscConstraintsIdx = _builder.add(ct);
                // padBloc(ct, col, blockZeroIdx);

                _builder.ensure_equality(miscConstraintsIdx, blockZeroIdx);
            }
    }

    // We are testing a bound where -d < x < d
    // we will need to test on positive numbers however, 0 < x + d < 2d => x + (D-d) < D
     std::vector<size_t> addBoundingConstraint(size_t degree, std::vector<FieldT> a, FieldT C, FieldT d, size_t log2D, size_t blocOneIdx) {
        std::vector<size_t> tmp;

        // Compute the number of blocs needed to process all degrees
        size_t numberOfBlocs = (degree/_blockSize);

        // For each bloc:
        for (size_t idx = 0; idx < numberOfBlocs; idx++) {

            std::vector<FieldT> vX, vXc(_blockSize, FieldT(0)), vY(_blockSize, FieldT(0));

            // Populate X and X + C for each bloc
            vX.assign(a.begin()+idx*_blockSize, a.begin()+(idx+1)*_blockSize);
            std::transform(vX.begin(), vX.end(), vY.begin(), [&](FieldT a){return a + C;});
            std::transform(vX.begin(), vX.end(), vXc.begin(), [&](FieldT a){return a + d;});

            // Bound Abs(X)
            size_t blocXidx = _builder.add(new variable_constraint<FieldT>());
            _variables.add(blocXidx,vX);
            tmp.push_back(blocXidx);

            size_t blocXAbsidx = _builder.add(new linear_constraint<FieldT>(
                    {
                        {blocXidx,blocOneIdx},
                        {FieldT(1),d}
                    }));

            addBitDecomposition(vXc, blocXAbsidx, log2D);

            // Bound Y = X + C
            size_t blocYidx = _builder.add(new linear_constraint<FieldT>(
                    {
                        {blocXidx, blocOneIdx},
                        {FieldT(1),C}
                    }));

            addBitDecomposition(vY, blocYidx, log2D);
        }

        return tmp;
    }

    // Logic for bit decomposition
    void addBitDecomposition(vector<FieldT> varVal, size_t varBlocIdx, size_t log2D) {

        // Linear constraints for bits
        std::vector<size_t> bitDecompositionBlocIdx;
        std::vector<FieldT> bitDecompositionScalar;

        // bits X0 .. XD
        for (size_t idx = 0; idx < log2D; idx++) {
            vector<FieldT> data;

            for (size_t col = 0; col < varVal.size(); col++) {
                uint64_t bit = static_cast<uint64_t>((static_cast<uint64_t>(varVal[col].getValue()) & (1ull << idx)) >> idx);

                // assert(bit<2);
                data.push_back(bit);
            }

            size_t bitBlocIdx = _builder.add(new variable_constraint<FieldT>());
            _variables.add(bitBlocIdx,data);

            // Quadratic constraints for bits
            // We don't need to create a variable for that, this constraint should equal 0
            size_t bitSquaredIdx = _builder.add(new quadratic_constraint<FieldT>(
                        {bitBlocIdx},
                        {bitBlocIdx})
                        );

            _builder.ensure_equality(bitBlocIdx, bitSquaredIdx);

            // Linear constraint
            bitDecompositionBlocIdx.push_back(bitBlocIdx);
            bitDecompositionScalar.push_back(FieldT(1ull << idx));
        }

        // Add the full number to the bit decomposition linear constraint
        size_t bitRecompositionIdx = _builder.add(new linear_constraint<FieldT>(
                    bitDecompositionBlocIdx,
                    bitDecompositionScalar)
                    );

        _builder.ensure_equality(bitRecompositionIdx, varBlocIdx);
    };


    std::vector<FieldT> fsPickEvaluationPoints() {
        size_t attempt = 0;

        size_t nbvalidEvaluationPts = 0;

        // Fiat-Shamir
        std::vector<FieldT> evaluationPoints;

        while (evaluationPoints.size() < _publicData.ptNumberEvaluation) {
            uint64_t rand; 

            // First we need to ensure that this evaluation point is not an (n+1) root of unity,
            // i.e. w^n = - 1, w^n + 1 = 0 
            if (this->_csPickRandom) {
                // Fiat-Shamir Randomness Pick
                rand = (*this->_randomness) & ((1ull << 62) -1); // The primes are 62-bit numbers so we can discard the last two bits
                this->_randomness++;
            } else {
                rand = 0;
            }

            // Testing that randomness is within the field to uphold uniform distribution
            if (rand < params<uint64_t>::P[this->_publicData.modulusIdx]) {
                FieldT evaluationPoint(rand);
                
                // Is this a valid evaluation point?
                if ((evaluationPoint^(_publicData.degree)) + FieldT(1) != FieldT(0)) {
                    evaluationPoints.emplace_back(evaluationPoint);
                }
            }

            attempt++;
            if (attempt >= _publicData.ptNumberEvaluation * evaluationPointsOverSamplingFactor) {
                throw internalError("exhausted randomness for picking evaluation points");
            }
        }

        return evaluationPoints;
    }

    std::pair<size_t, size_t> createBlocs0and1() {
        return {_builder.add(new constant_constraint<FieldT>(std::vector<FieldT>(_blockSize, FieldT(0)))), _builder.add(new constant_constraint<FieldT>(std::vector<FieldT>(_blockSize, FieldT(1))))};
    }

    public:
        expressNPStatement(PublicData& pdata, SecretData& sdata, builder<FieldT>& builder, vector<vector<FieldT>>* v, size_t blockSizeIn, uint64_t *&randomness) : _publicData(pdata), _secretData(sdata), _blockSize(blockSizeIn), _variables(v, blockSizeIn), _builder(builder), _randomness(randomness) {};

        // Build a Basic Constraint Set (Matrix Multiplication)
        void buildConstraintSet(Statements expressNPStatement, bool csRandomness = true) {

            this->_csPickRandom = csRandomness;

            switch(expressNPStatement) {
                case Statements::Generic_Addition:
                    {
                        // Reserve _variables input for all blocs
                        size_t row1idx;
                        size_t row2idx;

                        row1idx = _builder.add(new variable_constraint<FieldT>());
                        row2idx = _builder.add(new variable_constraint<FieldT>());

                        std::vector<FieldT> row1data = ligero::math::sampleRandomVectorOverField(this->_blockSize, FieldT(FieldT::getModulus()-1), static_cast<double>(584735));
                        std::vector<FieldT> row2data = ligero::math::sampleRandomVectorOverField(this->_blockSize, FieldT(FieldT::getModulus()-1), static_cast<double>(584735));

                        _variables.add(row1idx,row1data);
                        _variables.add(row2idx,row2data);

                        // Perform the addition
                        size_t rowOutIdx = _builder.add(new linear_constraint<FieldT>({
                                    row1idx,
                                    row2idx
                                }));
                    }
                    break;

                case Statements::Generic_Multiplication:
                    {
                        // Multiplying two random blocs
                        // Reserve _variables input for all blocs
                        size_t row1idx;
                        size_t row2idx;

                        row1idx = _builder.add(new variable_constraint<FieldT>());
                        row2idx = _builder.add(new variable_constraint<FieldT>());

                        std::vector<FieldT> row1data = ligero::math::sampleRandomVectorOverField(this->_blockSize, FieldT(FieldT::getModulus()-1), static_cast<double>(584735));
                        std::vector<FieldT> row2data = ligero::math::sampleRandomVectorOverField(this->_blockSize, FieldT(FieldT::getModulus()-1), static_cast<double>(584735));

                        _variables.add(row1idx,row1data);
                        _variables.add(row2idx,row2data);

                        // Perform the multiplication
                        size_t rowOutIdx = _builder.add(new quadratic_constraint<FieldT>(
                                    {row1idx},
                                    {row2idx}
                                    ));
                    }
                    break;

                case Statements::Generic_IntrablocBoundingVariables:
                    {
                        // Constant blocs
                        auto [blocZeroIdx, blocOneIdx] = createBlocs0and1();
                        
                        // Inputs
                        std::vector<FieldT> vars    = { 0, 2, 1, 7, 4, 5, 2, 1, 3};
                        std::vector<FieldT> ds      = { 8, 8, 8, 8, 8, 8, 8, 8, 8};

                        // Actual Constraint
                        addIntrablocBoundingConstraint(vars, ds, blocOneIdx, blocZeroIdx);

                    }
                    break;

                case Statements::Generic_InterblocBoundingVariables:
                    {
                        auto [d, C, D, log2D] = rephraseAs2Nct((1ull << 14)+(1ull << 2));

                        std::vector<FieldT> vX, vY;

                        for (size_t i = 0; i < _blockSize; i++) {
                            FieldT X((1ull << 11) + (1ull << 7) + (i % 5475));
                            vX.push_back(X);
                            vY.push_back(X + C);
                        }

                        // populating the constraints
                        // Constant Bloc C
                        size_t blocCidx = _builder.add(new constant_constraint<FieldT>(vector<FieldT>(_blockSize, C)));

                        // X
                        size_t blocXidx = _builder.add(new variable_constraint<FieldT>());
                        _variables.add(blocXidx,vX);
                        addBitDecomposition(vX, blocXidx, log2D);

                        // Y
                        DBG("X = " << vX[0].getValue());
                        DBG("Y = " << vY[0].getValue());
                        DBG("X + C = Y");

                        size_t blocYidx = _builder.add(new linear_constraint<FieldT>({
                                    blocXidx,
                                    blocCidx,
                                }));

                        addBitDecomposition(vY, blocYidx, log2D);
                    }
                    break;

                case Statements::Generic_IntraBloc:
                    {
                        std::vector<FieldT> a = ligero::math::sampleRandomVectorOverField(this->_blockSize, FieldT(FieldT::getModulus()-1), static_cast<double>(584735));
                        size_t blocAidx = _builder.add(new variable_constraint<FieldT>());
 
                        _variables.add(blocAidx,a);

                        // Add intrabloc constraint
                        transformation_constraint<FieldT> *newTransformationConstraint = new transformation_constraint<FieldT>();

                        for (size_t idx = 0; idx < 1+ 0*this->_blockSize/2; idx++) {
                            newTransformationConstraint->block.emplace_back(blocAidx);
                            newTransformationConstraint->scalar.emplace_back(FieldT(2));
                            newTransformationConstraint->source_position.emplace_back(idx);
                            newTransformationConstraint->target_position.emplace_back(idx + this->_blockSize/2);
                        }

                        _builder.add(newTransformationConstraint);
                    }
                    break;

                case Statements::Generic_PolyEval:
                    {
                        // Constant blocs
                        auto [blocZeroIdx, blocOneIdx] = createBlocs0and1();
                        size_t degree = 1ull << 16;
                        assert(_publicData.degree == degree);
                        
                        //_publicData.degree = degree;
                        
                        // Poly
                        // f(x) = a0
                        {
                            // inputs
                            std::vector<FieldT> basicPolyCoefs(degree,FieldT(0));
                            basicPolyCoefs[0] = FieldT(7);
                            vector<FieldT> evals = {4,7,35};

                            // expected values from the evaluation of the poly
                            size_t valuationIdx = (_builder.add(new variable_constraint<FieldT>()));
                            std::vector values(this->_blockSize, FieldT(0));

                            values[0] = basicPolyCoefs[0];
                            values[1] = basicPolyCoefs[0];
                            values[2] = basicPolyCoefs[0];
                            
                            _variables.add(valuationIdx, values);

                            vector<size_t> basicPolyIdxs = addSplitVariable(basicPolyCoefs);

                            // Is evaluation correct?
                            size_t blocCol = 0;
                            transformation_constraint<FieldT> *transformationConstraint = new transformation_constraint<FieldT>;

                            for (auto eval : evals) {
                                std::vector<FieldT> powersOfEval(degree);
                                for (size_t powIdx = 0; powIdx < degree; powIdx++) {
                                    powersOfEval[powIdx] = eval^powIdx;
                                }

                                addEvalPt(transformationConstraint, blocCol++, FieldT(1), basicPolyIdxs, powersOfEval);
                            }

                            // padBloc(transformationConstraint, blocCol, blocZeroIdx);
                            size_t evaluationBasicPolyIdx = _builder.add(transformationConstraint);

                            _builder.ensure_equality(evaluationBasicPolyIdx, valuationIdx);
                        }

                        // f(x) = a1 * x
                        {
                            // inputs
                            std::vector<FieldT> basicPolyCoefs(degree,FieldT(0));
                            basicPolyCoefs[1] = FieldT(1);
                            vector<FieldT> evals = {4,7,35};
                            // vector<size_t> res = {5,6,34};

                            // expected values from the evaluation of the poly
                            size_t valuationIdx = (_builder.add(new variable_constraint<FieldT>()));
                            std::vector values(this->_blockSize, FieldT(0));

                            values[0] = evals[0];
                            values[1] = evals[1];
                            values[2] = evals[2];
                            
                            _variables.add(valuationIdx, values);

                            vector<size_t> basicPolyIdxs = addSplitVariable(basicPolyCoefs);

                            // Is evaluation correct?
                            size_t blocCol = 0;
                            transformation_constraint<FieldT> *transformationConstraint = new transformation_constraint<FieldT>;

                            for (auto eval : evals) {
                                std::vector<FieldT> powersOfEval(degree);

                                for (size_t powIdx = 0; powIdx < degree; powIdx++) {
                                    powersOfEval[powIdx] = eval^powIdx;
                                }

                                addEvalPt(transformationConstraint, blocCol++, FieldT(1), basicPolyIdxs, powersOfEval);
                            }

                            // padBloc(transformationConstraint, blocCol, blocZeroIdx);
                            size_t evaluationBasicPolyIdx = _builder.add(transformationConstraint);

                            _builder.ensure_equality(evaluationBasicPolyIdx, valuationIdx);
                        }
                    }
                    break;


                case Statements::RSARound3_keyGen:
                    {
                        // Preprocessing: Derive & store phis
                        Preprocessing();
                        
                        // Constant blocs
                        auto [blocZeroIdx, blocOneIdx] = createBlocs0and1();
                        
                        // NP portions
                        NP_RSARound3_keyGen(blocZeroIdx, blocOneIdx);

                        // Clean up
                        free(this->_phis);
                    }
                    break;

                case Statements::RSARound4_preSieving:
                    {
                        // Preprocessing: Derive & store phis
                        Preprocessing();

                        // Constant blocs
                        auto [blocZeroIdx, blocOneIdx] = createBlocs0and1();

                        // NP portions
                        NP_RSARound4_preSieving(blocZeroIdx, blocOneIdx);

                        // Clean up
                        free(this->_phis);
                    }
                    break;

                case Statements::RSARound5_preSieving:
                    {

                        // Constant blocs
                        auto [blocZeroIdx, blocOneIdx] = createBlocs0and1();

                        // NP portions
                        NP_RSARound5_preSieving(blocZeroIdx, blocOneIdx);

                        // Clean up
                        free(this->_phis);
                    }
                    break;

                case Statements::RSARound6_partialDecrypt:
                    {
                        // Preprocessing: Derive & store phis
                        Preprocessing();

                        // Constant blocs
                        auto [blocZeroIdx, blocOneIdx] = createBlocs0and1();

                        // NP portions
                        NP_RSARound6_partialDecrypt(blocZeroIdx, blocOneIdx);

                        // Clean up
                        free(this->_phis);
                    }
                    break;

                case Statements::RSARound7_candidateGeneration:
                    {
                        // Preprocessing: Derive & store phis
                        Preprocessing();

                        // Constant blocs
                        auto [blocZeroIdx, blocOneIdx] = createBlocs0and1();

                        // NP portions
                        NP_RSARound7_candidateGeneration(blocZeroIdx, blocOneIdx);

                        // Clean up
                        free(this->_phis);
                    }
                    break;

                case Statements::RSARound8_beaversTriples:
                    {
                        // Preprocessing: Derive & store phis
                        Preprocessing();

                        // Constant blocs
                        auto [blocZeroIdx, blocOneIdx] = createBlocs0and1();

                        // NP portions
                        NP_RSARound8_beaversTriples(blocZeroIdx, blocOneIdx);

                        // Clean up
                        free(this->_phis);
                    }
                    break;

                case Statements::RSARound10_jacobiTest:
                    {
                        // v = gamma, u = gamma ^ x mod N
                        // x = log_v(u) mod N
                        //
                        // Public data:
                        // v = gamma - public 
                        //
                        // Secret data:
                        // x - secret
                        // r - secret
                        // a - secret a = v ^ r
                        //
                        // in this step:
                        // here goes Jacobi test verification
                        // v = gamma, u = gamma ^ x mod N
                        // x = log_v(u) mod N
                        // a = v ^ r mod N
                        // 1. Pick random r, a = v ^ r mod N, send a
                        // 2. Verifier sends message e
                        // 3. Prover sends z = r + e * x (NP statement)
                        // 4. Verifier checks v ^ z = a * (u ^ e) (NOT NP statement)
                        //
                        // e' = hash(a, v, u) 
                        // e = prg(seed=e')
                        // 3. Proof: a, e, z -> v^z = a * u ^ e 

                            // Constant blocs
                        auto [blocZeroIdx, blocOneIdx] = createBlocs0and1();
                        NP_RSARound10_jacobiTest(blocZeroIdx, blocOneIdx);

                        // addEquationJacobi( 
                        //     transformationConstraint,
                        //     evaluationPoints,
                        //     blocZeroIdx,
                        //     blocOneIdx,
                        //     1,// _publicData.gamma,
                        //     blocCol);

                    }
                    break;

                case Statements::RSARound11_12_jacobiGCD:
                    {
                        // Preprocessing: Derive & store phis
                        Preprocessing();

                        // Constant blocs
                        auto [blocZeroIdx, blocOneIdx] = createBlocs0and1();

                        // NP portions
                        NP_RSARound11_12_jacobiGCD(blocZeroIdx, blocOneIdx);

                        // Clean up
                        free(this->_phis);
                    }
                    break;

                case Statements::RSACeremony_Connecting_Proofs:
                    {
                        // Pre-requisites: Steps 4, 5, 7, 8 and 11/12

                        // Preprocessing: Derive & store phis
                        Preprocessing();

                        // Constant blocs
                        auto [blocZeroIdx, blocOneIdx] = createBlocs0and1();

                        // NP portions
                        NP_RSARound4_preSieving(blocZeroIdx, blocOneIdx);
                        NP_RSARound5_preSieving(blocZeroIdx, blocOneIdx);
                        NP_RSARound7_candidateGeneration(blocZeroIdx, blocOneIdx);
                        NP_RSARound8_beaversTriples(blocZeroIdx, blocOneIdx);
                        NP_RSARound11_12_jacobiGCD(blocZeroIdx, blocOneIdx);
                        NP_Connecting_Proofs(blocZeroIdx, blocOneIdx);

                        // Clean up
                        free(this->_phis);
                    }
                    break;

                case Statements::RSACeremony_Equating_Variables:
                    {
                        // Pre-requisites: Steps 7, 8 and 11/12

                        // Preprocessing: Derive & store phis
                        Preprocessing();

                        // Constant blocs
                        auto [blocZeroIdx, blocOneIdx] = createBlocs0and1();

                        // NP portions
                        NP_RSARound7_candidateGeneration(blocZeroIdx, blocOneIdx);
                        NP_RSARound8_beaversTriples(blocZeroIdx, blocOneIdx);
                        NP_RSARound11_12_jacobiGCD(blocZeroIdx, blocOneIdx);
                        NP_Equate_Variables(blocZeroIdx, blocOneIdx);

                        // Clean up
                        free(this->_phis);
                    }
                    break;

                case Statements::RSACeremony_Bounding_Variables:
                    {
                        // Preprocessing: Derive & store phis
                        Preprocessing();

                        // Constant blocs
                        auto [blocZeroIdx, blocOneIdx] = createBlocs0and1();

                        // NP portions
                        NP_RSACeremony_Bounding_Variables(blocZeroIdx, blocOneIdx);

                        // Clean-up
                        free(this->_phis);
                    }
                    break;

                case Statements::RSACeremony:
                    {
                        // Preprocessing: Derive & store phis
                        Preprocessing();

                        // Constant blocs
                        auto [blocZeroIdx, blocOneIdx] = createBlocs0and1();

                        NP_RSARound3_keyGen(blocZeroIdx, blocOneIdx);
                        NP_RSARound4_preSieving(blocZeroIdx, blocOneIdx);
                        NP_RSARound5_preSieving(blocZeroIdx, blocOneIdx);
                        NP_RSARound6_partialDecrypt(blocZeroIdx, blocOneIdx);
                        NP_RSARound7_candidateGeneration(blocZeroIdx, blocOneIdx);
                        NP_RSARound8_beaversTriples(blocZeroIdx, blocOneIdx);
                        NP_RSARound11_12_jacobiGCD(blocZeroIdx, blocOneIdx);
                        NP_Connecting_Proofs(blocZeroIdx, blocOneIdx);
                        NP_Equate_Variables(blocZeroIdx, blocOneIdx);

                        // Clean-up
                        free(this->_phis);
                    }
                    break;

                case Statements::RSACeremony_Rounds3to6:
                    {
                        // Preprocessing: Derive & store phis
                        Preprocessing();

                        // Constant blocs
                        auto [blocZeroIdx, blocOneIdx] = createBlocs0and1();

                        NP_RSARound3_keyGen(blocZeroIdx, blocOneIdx);
                        NP_RSARound4_preSieving(blocZeroIdx, blocOneIdx);
                        NP_RSARound5_preSieving(blocZeroIdx, blocOneIdx);
                        NP_RSARound6_partialDecrypt(blocZeroIdx, blocOneIdx);

                        // Clean-up
                        free(this->_phis);
                    }
                    break;

                case Statements::RSACeremony_EarlyCommitments:
                    {
                        std::vector<FieldT> ei;
                        std::vector<FieldT> si;
                        std::vector<FieldT> xi;

                        // Secret Data: either source or initialize if not available
                        if (_secretData.isAvailable) {
                            ei = embed(_secretData.eiP, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
                            si = embed(_secretData.siP, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
                            xi = embed(_secretData.xi, _publicData.modulusIdx * _publicData.degree, _publicData.degree);
                        } else {
                            ei.assign(_publicData.degree, FieldT(0));
                            si.assign(_publicData.degree, FieldT(0));
                            xi.assign(_publicData.degree, FieldT(0));
                        }

                        this->_siCoefsBlocIds = addSplitVariable(si);
                        this->_eiCoefsBlocIds = addSplitVariable(ei);
                        this->_xiCoefsBlocIds = addSplitVariable(xi);
                    }
                    break;

                case Statements::GenerateData:
                default:
                    break;
            }
        }
    };
