# ISOCC'23-Fault bound

# Author

**Wonyeong Jung** 
- Email: jwy6617@gmail.com

**Dongwhee Kim**
- Email: xyz12347976@gmail.com
- Google Scholar: https://scholar.google.com/citations?user=8xzqA8YAAAAJ&hl=ko&oi=ao

Paper title: [ISOCC'23] Fault_bound

Paper URL: **(to be updated later)**

# Overview
![An Overview of the Fault_bound](https://github.com/xyz123479/ISOCC_23-Fault-Bound/blob/main/Fault_bound.png)

# Code flows (Fault_sim.cpp)
- 1. Reading OD-ECC, RL-ECC H-Matrix.txt: It's fine not to use RL-ECC H-Matrix.txt.
- 2. Setting output function name: output.S file.
- 3. **(Start loop)** DDR5 ECC-DIMM setup
- 4. Initialize all data in 10 chips to 0: Each chip has 136 bits of data + redundancy.
- 5. Error injection: Errors occur based on the following probabilities:
>> SE: 40%, DE: 30%, SCE: 14%, SE+SE: 16%
- 6. **(Fill in the code)** Apply OD-ECC: Implementation
>> Apply the Hamming SEC code of (136, 128) to each chip.

>> After running OD-ECC, the redundancy of OD-ECC does not come out of the chip (128bit data).
- 7. **(Fill in the code)** Apply RL-ECC
>> Run (80, 64) RL-ECC by bundling two beats.
>> Please feel free to use any ECC code.
>> 16 Burst Length (BL) creates one memory transfer block (64B cacheline + 16B redundancy).
>> In DDR5 x4 DRAM, because of internal prefetching, only 64bit of data from each chip's 128bit data is actually transferred to the cache.
>> For this, create two memory transfer blocks for 128-bit data and compare them.
- 8. Report CE/DUE/SDC results.
- 9. **(End loop)** Derive final results.

# DIMM configuration (per-sub channel)
- DDR5 ECC-DIMM
- Num of rank: 1
- Beat length: 40 bit
- Burst length: 16
- Num of data chips: 8
- Num of parity chips: 2
- Num of DQ: 4 (x4 chip)

# ECC configuration
- OD-ECC: (136, 128) Hamming SEC code **[1]** 'or' SEC code with bounded_Fault **[2]**
- RL-ECC: Chipkill-correct ECC **[3-4]**

# Error pattern configuration
- SE(SBE): per-chip Single Bit Error
- DE(DBE): per-chip Double Bit Error
- CHIPKILL(SCE): Single Chip Error (All Random)

# Error Scenario configuration
- SE(SBE): Among 10 chips, there's a single bit error (SE[Single Bit Error]) occurring in just one chip, with the remaining 9 chips having no errors
- DE(DBE): Among 10 chips, there's a double bit error (DE[Double Bit Error]) occurring in just one chip, with the remaining 9 chips having no errors
- CHIPKILL(SCE): Among 10 chips, there's a random error (SCE [Single Chip Error]) occurring in just one chip, with the remaining 9 chips having no errors. Errors can occur up to a maximum of 136 bits
- SE(SBE)+SE(SBE): Among 10 chips, there's a single bit error (SE[Single Bit Error]) occurring in each of two chips, with the remaining 8 chips having no errors

# References
- **[1]** Hamming, Richard W. "Error detecting and error correcting codes." The Bell system technical journal 29.2 (1950): 147-160.
- **[2]** Criss, Kjersten, et al. "Improving memory reliability by bounding DRAM faults: DDR5 improved reliability features." The International Symposium on Memory Systems. 2020.
- **[3]** Reed, Irving S., and Gustave Solomon. "Polynomial codes over certain finite fields." Journal of the society for industrial and applied mathematics 8.2 (1960): 300-304.
- **[4]** https://www.amd.com/system/files/TechDocs/42301_15h_Mod_00h-0Fh_BKDG.pdf
- **[5]** Pontarelli, Salvatore, et al. "Low delay single symbol error correction codes based on reed solomon codes." IEEE transactions on computers 64.5 (2014): 1497-1501.

