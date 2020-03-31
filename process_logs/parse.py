#!/usr/local/bin/python3.7

from os import path
from datetime import datetime
import re
import glob
from statistics import mean, stdev

### REGULAR EXPRESSIONS FOR PARSING PARTS OF LOGS

logts = re.compile('([0-9]{2}):([0-9]{2}):([0-9]{2}).([0-9]{3})')

# 21:12:23.981 ‹ ,1.a. Overall speed, , ,374.32 MB,374.32 MB,06:21.692612,06:21.692616,
s1a = re.compile('1\.a\. Overall speed, , ,([0-9\.]+) MB,([0-9\.]+) MB,(\d\d):(\d\d)\.(\d+)') #, , ,([0-9\.]+) MB,([0-9\.]+) MD,(\d\d):(\d\d)\.(\d+)')

# 21:06:02.288 Registration completed for 2 out of 2
reg = re.compile('Registration completed')
preg = re.compile('registering with coordinator')
# 21:06:24.001 ‹RSA Ceremony, , , ,304.11 MB,304.11 MB,00:21.713310,00:21.713311,
# 21:25:11.714 Found 1 valid moduli:
# 07:41:24.961 No candidates found.
passive = re.compile('(Found . valid moduli)|(No candidates found)')

# 21:12:23.981 Verified all proofs; ceremony successful.
allver = re.compile('Verified all proofs')

# 21:06:03.498 MessageType: PUBLIC_KEY_A_VALUE message size: 11010105 bytes
msgtype = re.compile('MessageType: ([A-Z\_]+) message size: ([0-9]+) bytes')

memory = re.compile('Peak Memory = ([0-9]+) Kb')

# 21:26:08.706 Verifying modulus idx:0
vstart = re.compile('Verifying modulus idx:([0-9]+)')
vend   = re.compile('Verification for modulus idx:[0-9]+ succeeded')

sendproof = re.compile('Send Proof for Modulus ([0-9]+)')

# MSGS should occur in this order. Todo: check for this in each log
msg_list = [
  'PROTOCOL_CONFIG',
  'PUBLIC_KEY_A_VALUE',
  'PUBLIC_KEY_B_VALUE',
  'ASSIGNMENT_PN',
  'ENCRYPTED_X_VALUE',
  'ENCRYPTED_XY_PLUS_Z_VALUE',
  'PS_SIEVING_FLAGS',
  'AX_BY_VALUE',
  'MODULUS_CANDIDATE',
  'POST_SIEVE',
  'GAMMA_SHARES',
  'GAMMA_RANDOM_SEED_VALUE',
  'DISCARD_FLAGS',
  'GCD_RAND_SHARES',
  'AX_BY_VALUE',
  'DISCARD_FLAGS',
  'FOUND_MODULI',
]

active_list = [
  'GATHER_PUBLIC_DATA',
  'GATHER_PROOF_0',
  'GATHER_PROOF_1',
  'GATHER_PROOF_2',
  'GATHER_PROOF_3',
  'GATHER_PROOF_4',
  'GATHER_PROOF_5',
  'GATHER_PROOF_6',
  'GATHER_PROOF_7',
  'GATHER_PROOF_8',
  'GATHER_PROOF_9',
  'GATHER_PROOF_10',
  'GATHER_PROOF_11',
  'GATHER_PROOF_12',
  'GATHER_PROOF_13',
  'GATHER_PROOF_14',
  'GATHER_PROOF_15',
  'GATHER_PROOF_16',
  'GATHER_PROOF_17',
  'GATHER_PROOF_18',
  'GATHER_PROOF_19',
  'GATHER_PROOF_20',
]

party_msg_list = [
#'ID_PARTY',
'PUBLIC_KEY_A_SHARES',
'PUBLIC_KEY_B_SHARES',
'ENCRYPTED_X_SHARES',
'ENCRYPTED_XY_PLUS_Z_SHARES',
'PARTIAL_XY_MINUS_Z_SHARES',
'AX_BY_SHARES',
'AXB_MINUS_BYA_SHARES',
'MUTHU_ACK',
'GAMMA_RANDOM_SEED_SHARES',
'EXPONENTIATED_GAMMA_VALUE',
'GCD_AX_BY_SHARES',
'AXB_MINUS_BYA_SHARES'
]

class Experiment:
    def __init__(self, name):
      self.name = name
      self.registration = []
      self.party_registration = []
      self.passive = []
      self.active = []
      self.overall = []
      self.msg_ts = {}
      self.msg_sz = {}
      self.party_msg_ts = {}
      self.party_msg_sz = {}
      self.memory = []
      self.vidle = {}   # verifier idle times
      self.vwork = {}   # verifier work times

    def summary(self):
      print(self.name + '  ' + str(len(self.registration)) + ' runs')
      print(f'  registration: {mean(self.registration):.2f}')
      print(f'p registration: {mean(self.party_registration):.2f}')
      print(f'       passive: {self.avg_passive():.2f}')
      print(f'        active: {self.avg_active():.2f}')
      for m in msg_list:
        print(f'  {m.rjust(26)} C: {self.avg_msg(m):.2f}')
      for m in party_msg_list:
        print(f'  {m.rjust(26)} P: {self.avg_party_msg(m):.2f}')

      # for m in active_list:
      #   print(f'  {m.rjust(26)} C: {mean(self.msg_ts[m]):.2f}')

      # for m in active_list[1:]:
      #   print(f'  {m.rjust(26)} P: {mean(self.party_msg_ts[m]):.2f}')

    def avg_passive(self):
      return mean(self.passive)
    def std_passive(self):
      if len(self.passive)>1:
        return stdev(self.passive)
      else:
        return 0


    def avg_active(self):
      if len(self.active) > 0:
        return mean(self.active)
      else:
        return -1
    def std_active(self):
      if len(self.active) > 1:
        return stdev(self.active)
      else:
        return -1


    def avg_reg(self):
      return mean(self.registration)
    def avg_mem(self):
      if len(self.memory) > 0:
        return mean(self.memory)
      else:
        return -1
    def avg_msg(self, msg):
      if msg in self.msg_ts.keys() and len(self.msg_ts[msg]) > 0:
        return mean(self.msg_ts[msg])
      else:
        return -1
    def avg_party_msg(self, msg):
      if msg in self.party_msg_ts.keys() and len(self.party_msg_ts[msg]) > 0:
        return mean(self.party_msg_ts[msg])
      else:
        return -1
    def avg_msg_sz(self, msg):
      if msg in self.msg_sz.keys() and len(self.msg_sz[msg]) > 0:
        return mean(self.msg_sz[msg])
      else:
        return -1

    def std_msg_sz(self, msg):
      if msg in self.msg_sz.keys() and len(self.msg_sz[msg]) > 0:
        return stdev(self.msg_sz[msg])
      else:
        return 0

    def avg_party_msg_sz(self, msg):
      if msg in self.party_msg_sz.keys() and len(self.party_msg_sz[msg]) > 0:
        return mean(self.party_msg_sz[msg])
      else:
        return -1

    def std_party_msg_sz(self, msg):
      if msg in self.party_msg_sz.keys() and len(self.party_msg_sz[msg]) > 0:
        return stdev(self.party_msg_sz[msg])
      else:
        return 0

    def std_party_msg(self, msg):
      if msg in self.party_msg_ts.keys() and len(self.party_msg_ts[msg]) > 0:
        return stdev(self.party_msg_ts[msg])
      else:
        return -1

    def avg_vidle(self, idx):
      if len(self.vidle) > 0:
        return mean(self.vidle[idx])
      else:
        return -1

    def avg_vwork(self, idx):
      if len(self.vwork) > 0:
        return mean(self.vwork[idx])
      else:
        return -1

    def std_vidle(self, idx):
      if len(self.vidle) > 0:
        return stdev(self.vidle[idx])
      else:
        return -1

    def std_vwork(self, idx):
      if len(self.vwork) > 0:
        return stdev(self.vwork[idx])
      else:
        return -1



def ts(line):
  m = logts.match(line)
  ts = 0
  if m:
    ts = int(m.group(1))*60*60*1000 + int(m.group(2))*60*1000 + int(m.group(3))*1000 + int(m.group(4))
  return ts


def coordinator_parser(rp, exp):
  filepath = rp + '/coordinator.log'
  cnt = 0
  with open(filepath) as fp:
    line = fp.readline()
    start = ts(line)
    prot_start = start
    last_sent_msg = start
    gp = 0   
    tt = 0
    idx = 0
    gp_time = 0
    while line:
      #print("Line {}: {}".format(cnt, line.strip()))
      t = ts(line)
      m = reg.search(line)
      if m:
        print('   ' + str(t-start) + ' registration')
        exp.registration.append(t-start)
        prot_start = t
        last_sent_msg = prot_start
        print(f'   setting prot_start {prot_start}')

      m = passive.search(line)
      if m:
        print(f'   {(t-prot_start)}  passive done. sum:{tt}')
        exp.passive.append(t-prot_start)
        gp_start = t

      m = allver.search(line)
      if m:
        print('   ' + str(t-prot_start) + ' all done')
        exp.active.append(t-prot_start)

      m = msgtype.search(line)
      if m:
        key = m.group(1)
        if key == "GATHER_PROOFS":
          gp = gp + 1
          key = "GATHER_PROOF_" + str(gp)
#        print('    ' + str(t-last_sent_msg) + ' ' + key + ' ' + m.group(2))
        exp.msg_sz.setdefault(key, []).append(int(m.group(2)))

        # handle proof timing below
        if m.group(1) != "GATHER_PROOFS":
          exp.msg_ts.setdefault(key, []).append(t-last_sent_msg)

        tt += (t-last_sent_msg)
        last_sent_msg = t

      # gather proofs are sent in batches of n to verifiers in order.
      m = sendproof.search(line)
      if m:
        nidx = int(m.group(1))
        if nidx > idx:    # moving on to next proof
          key = "GATHER_PROOF_" + str(idx)
          exp.msg_ts.setdefault(key, []).append(t-gp_start)
          gp_start = t
          idx = nidx

      m = s1a.search(line)
      if m:
        dur = int(m.group(3))*60*1000 + int(m.group(4))*1000 + int(m.group(5))/1000
 #       print('   ' + str(t-prot_start) + ' Overall: ' + m.group(1) + ' ' + m.group(2) + ' ' + str(dur))
        exp.overall.append(dur)


      #print(str(t) + "\n")
      # m = s1a.match(line)
      # if m:
      #   print(line)
      #   print(m.group(1) + ' ' + m.group(2) )
      line = fp.readline()
      cnt += 1


def party_parser(rp, exp):
  cnt = 0
  for f in glob.iglob(rp+'/party_full_protocol_*.log'):
    cnt = cnt + 1
#    print('    ' + f)
    with open(f) as fp:
      line = fp.readline()
      start = ts(line)
      last_sent_msg = start
      gp = 0    # the GATHER_PROOF counter, since we want to separate by party  
      while line:
        t = ts(line)

        # reset start as soon as reg is done
        m = preg.search(line)
        if m:
          exp.party_registration.append(t-start)
          start = t
          last_sent_msg = t

        ## parsing msgtype so far
        m = msgtype.search(line)
        if m:
          key = m.group(1)
          if key == "GATHER_PROOFS":
            key = "GATHER_PROOF_" + str(gp)
            gp = gp + 1
#          print('    ' + str(t-last_sent_msg) + ' ' + key + ' ' + m.group(2))
          exp.party_msg_ts.setdefault(key, []).append(t-last_sent_msg)
          exp.party_msg_sz.setdefault(key, []).append(int(m.group(2)))
          last_sent_msg = t

        m = memory.search(line)
        if m:
          exp.memory.append(int(m.group(1)))

        line = fp.readline()

  print(f'   parsed {cnt} party log files')


def verifier_parser(rp, exp):
  for f in glob.iglob(rp+'/distributed_verifier*.log'):
    with open(f) as fp:
      line = fp.readline()
      last_start = ts(line)
      last_end = last_start
      idx = 0     # idx we are trying to verify next
      while line:
        t = ts(line)

        m = vstart.search(line)
        if m:
          last_start = ts(line)
          if idx > 0:
            idle = t-last_end
            # print(f'   idle {idle} {line}')
            exp.vidle.setdefault(idx,[]).append(idle)

        m = vend.search(line)
        if m:
          vtime = t-last_start
          exp.vwork.setdefault(idx,[]).append(vtime)
          idx = idx+1
          last_end = t
          # print(f'  work {vtime} {line}')

        line = fp.readline()


dirs = [ "2", "5", "10", "20", "50", "100", ]#"200", "500", "1000", "2000", "4046"  ]
maindir = './data/03-04-20/'

experiments = {}

for d in dirs:
    cnt = 1
    exp = Experiment(d)
    experiments[d] = exp
    while True:
        rp = maindir + d + '/run' + str(cnt)
        if path.exists(rp):
            print('Parsing %s' % rp)
            coordinator_parser(rp, exp)
            party_parser(rp, exp)
            verifier_parser(rp, exp)
        else:
            break
        cnt = cnt + 1
    print('   ['+d+'] Done parsing ' + str(cnt-1) + ' runs')
    exp.summary()


# make the summary tables of passive/active total time
print(f'#n   passive  active  reg')
for d in dirs:
  e = experiments[d]
  print(f'{d.ljust(5)} & {e.avg_passive()/1000:.1f} {e.std_passive()/1000:.1f} & {e.avg_active()/1000:.1f} {e.std_active()/1000:.1f}& & {e.avg_reg()/1000:.1f}\\\\[2pt] % {len(e.passive)} {len(e.active)} ')

print()
# make the per-message table
print(f'n   ', end=' ') 
for m in msg_list:
  print(f'{m}', end=' ')
print('PASS SUM')
for d in dirs:
  e = experiments[d]
  print(f'{d.ljust(5)}', end =" ")
  sum = 0
  for m in msg_list:
    sum += (e.avg_msg(m)/1000)
    print(f'{e.avg_msg(m)/1000:5.1f}', end =' ')
  print(f'{e.avg_passive()/1000:5.1f} {sum:5.1f}')


########## make the verifier idle/work time
print()
print(f'n  work0 idle1 work1 ... (for verifier.txt)') 
for d in dirs:
  e = experiments[d]
  print(f'{d.ljust(5)} {e.avg_vwork(0)/1000:.1f}', end=' ')
  for idx in range(1,21):
    print(f'{e.avg_vidle(idx)/1000:.1f} {e.avg_vwork(idx)/1000:.1f}', end=' ')
  print()

###### coordinator proof, for m2.txt
print()
for d in dirs:
  e = experiments[d]
  print(f'{d.ljust(5)}', end =" ")
  sum = 0
  for m in active_list[1:]:
    print(f'{e.avg_party_msg(m)/1000:5.2f}', end =' ')  
  print()


### msg size analysis

sz_list = [
'PUBLIC_KEY_A_VALUE',
'PUBLIC_KEY_B_VALUE',
'ENCRYPTED_X_VALUE',
'ENCRYPTED_XY_PLUS_Z_VALUE',
'PS_SIEVING_FLAGS',
'AX_BY_VALUE',
'MODULUS_CANDIDATE',
'AX_BY_VALUE',
]


print()
for m in sz_list:
  print(f'{m}', end=' ')
print('OTHER')
for d in dirs:
    e = experiments[d]
    print(f'{d.ljust(5)}', end =" ")
    other = e.avg_msg_sz('PROTOCOL_CONFIG') + e.avg_msg_sz('ASSIGNMENT_PN') + e.avg_msg_sz('POST_SIEVE') +e.avg_msg_sz('GAMMA_SHARES') + e.avg_msg_sz('GAMMA_RANDOM_SEED_VALUE') +e.avg_msg_sz('DISCARD_FLAGS') + e.avg_msg_sz('GCD_RAND_SHARES') +e.avg_msg_sz('DISCARD_FLAGS') + e.avg_msg_sz('FOUND_MODULI')
    for m in sz_list:
      print(f'{e.avg_msg_sz(m):5.1f}', end =' ')
    print(f'{other:5.1f}')

# echo "filename,idle,compute" >> $out
# for dataFile in $(ls data/distributed_*) 
# do
# idle=`cat $dataFile | grep '6.a. verification' | grep idle | tail -n +2  | awk -F',' '{ print $9 }' | python3 -c "import sys; x=sys.stdin.read().split('\n')[:-1]; print(sum([float(i.split(':')[0]) * 60 + float(i.split(':')[1]) for i in x]))"`
# compute=`cat $dataFile | grep '6.a. verification' | grep compute | awk -F',' '{ print $9 }' | python3 -c "import sys; x=sys.stdin.read().split('\n')[:-1]; print(sum([float(i.split(':')[0]) * 60 + float(i.split(':')[1]) for i in x]))"`
# echo "$dataFile,$idle,$compute" >> $out
# done
