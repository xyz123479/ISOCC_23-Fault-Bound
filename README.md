# ISOCC'23-Fault bound

# Author

**Wonyeong Jung** 
- Email: jwy6617@gmail.com

**Dongwhee Kim**
- Email: xyz12347976@gmail.com
- Google Scholar: https://scholar.google.com/citations?user=8xzqA8YAAAAJ&hl=ko&oi=ao

Paper title: [ISOCC'23] Synergistic Integration: An Optimal Combination of On-Die and Rank-Level ECC for Enhanced Reliability

Paper URL: **(to be updated later)**

# Overview
![An Overview of the Fault_bound](https://github.com/xyz123479/ISOCC_23-Fault-Bound/blob/main/Fault_bound.png)

# Code flows (Fault_sim.cpp)
- 1. Reading OD-ECC, RL-ECC H-Matrix.txt: It's fine not to use RL-ECC H-Matrix.txt.
- 2. Setting output function name: output.S file.
- 3. **(Start loop)** DDR5 ECC-DIMM setup
- 4. Initialize all data in 10 chips to 0: Each chip has 136 bits of data + redundancy.
- 5. Error injection: Scenario-based Error injection
- 6. Apply OD-ECC: Implementation
>> Apply the Hamming SEC code of (136, 128) to each chip.

>> After running OD-ECC, the redundancy of OD-ECC does not come out of the chip (128bit data).
- 7. Apply RL-ECC
>> 16 Burst Length (BL) creates one memory transfer block (64B cacheline + 16B redundancy).

>> In DDR5 x4 DRAM, because of internal prefetching, only 64bit of data from each chip's 128bit data is actually transferred to the cache.

>> For this, **create two memory transfer blocks for 128-bit data and compare them.**
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
- RL-ECC: Chipkill-correct ECC using RS (Reed-Solomon) code **[3]**
>> Applying **Restrained mode [4]**

# Error pattern configuration (2 chip errors)
- SE: per-chip Single bit Error
- MBBE: per-chip Muli bit Bounded Error
- SW: per-chip Single Word Error (all 4 bits)
- SP: per-chip Single Pin Error (More than 2 bits)
- CHIPKILL: Single Chip error (All random)

# Getting Started
- $ make clean
- $ make
- $ python run.py

# References
- **[1]** Hamming, Richard W. "Error detecting and error correcting codes." The Bell system technical journal 29.2 (1950): 147-160.
- **[2]** Criss, Kjersten, et al. "Improving memory reliability by bounding DRAM faults: DDR5 improved reliability features." The International Symposium on Memory Systems. 2020.
- **[3]** Reed, Irving S., and Gustave Solomon. "Polynomial codes over certain finite fields." Journal of the society for industrial and applied mathematics 8.2 (1960): 300-304.
- **[4]** Kim, Dongwhee, et al. "Unity ECC: Unified Memory Protection Against Bit and Chip Errors." Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis. 2023.

