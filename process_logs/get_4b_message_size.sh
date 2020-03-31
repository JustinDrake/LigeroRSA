out="raw/4.b.MessageSizePerTypeForCoordinator.csv"
echo "MessageType,Size(bytes)" >> $out
for MESSAGE_TYPE in "PROTOCOL_CONFIG" "PUBLIC_KEY_A_VALUE" "PUBLIC_KEY_B_VALUE" "ASSIGNMENT_P1" "ASSIGNMENT_PN" "ENCRYPTED_X_VALUE" "ENCRYPTED_XY_PLUS_Z_VALUE" "PS_SIEVING_FLAGS" "AX_BY_VALUE" "MODULUS_CANDIDATE" "POST_SIEVE" "GAMMA_SHARES" "GAMMA_RANDOM_SEED_VALUE" "DISCARD_FLAGS" "GCD_RAND_SHARES" "AX_BY_VALUE" "DISCARD_FLAGS"
do
cat data/coordinator.log | grep 'message size' | awk '{ print $3","$6 }' | grep ${MESSAGE_TYPE} >> $out
done
