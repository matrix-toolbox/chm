#!/usr/bin/bash
# 2022-07-30

# machine with 4GB of RAM can handle max. 40000 of CHUNK
# ~40 minutes on Intel(R) Core(TM) i3-5005U CPU @ 2.00GHz with Linux 0x00 5.18.11-arch1-1  

# machine with 4GB of RAM can handle max. 40000 of CHUNK 
# ~25 minutes on Intel(R) Core(TM) i5-2500 CPU @ 3.30GHz with Linux 0x05 5.16.2-arch1-1

# machine with 16GB of RAM can handle max. 100000 of CHUNK
# ~4 hours on Intel(R) Core(TM) i7-4770 CPU @ 3.40GHz  with micro$hit windows7

./BHrc.ELF SPLIT_FILEaa 12 40000
octave get_defect.m > defect_BH_12_12aa.txt
shred -zuv SPLIT_FILEaa

./BHrc.ELF SPLIT_FILEab 12 40000
octave get_defect.m > defect_BH_12_12ab.txt
shred -zuv SPLIT_FILEab

./BHrc.ELF SPLIT_FILEAac 12 40000
octave get_defect.m > defect_BH_12_12ac.txt
shred -zuv SPLIT_FILEac

./BHrc.ELF SPLIT_FILEad 12 40000
octave get_defect.m > defect_BH_12_12ad.txt
shred -zuv SPLIT_FILEad

# ...

