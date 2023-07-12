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
- 1. Setting ECC-block configuration & error scenarios.
- 2. Setting output function name: output.S file.
- 3. **(Start loop)** HBM2E ECC-block setup
- 4. Initialize all data in ECC-block to 0
- 5. Error injection: Errors occur based on the error scenarios. **(Caution!) This evaluation has no fault!**
- 6. Apply **OD-ECC (On-Die ECC)**
>> Prior work: Apply the Hsiao SEC-DED code of (104, 96) to each ECC block.
>> 
>> EPA-ECC: Apply the RS SSC-DSD code of [39, 36] to an ECC block.
- 7. Report CE/DUE/SDC results.
- 8. **(End loop)** Derive final results.

# DDR5 ECC-DIMM configuration [1]
- Data: 256 bit
- System ECC redundancy: 32 bit
- On-Die ECC redundancy: 24 bit
- Num of DQ: 64 (Psuedo-channel mode)
- Num of Redundancy-DQ: 8
- Burst Length: 4

# Getting Started
- $ make clean
- $ make
- $ python run.py

# Answer (.S files)
- CE: detected and corrected error
- DUE: detected but uncorrected error
- SDC: Silent Data Corruption

# RUN_NUM is in Fault_sim.cpp file (iteration count)
- #define RUN_NUM 100000000

# References
- **[1]** Chun, Ki Chul, et al. "A 16-GB 640-GB/s HBM2E DRAM with a data-bus window extension technique and a synergetic on-die ECC scheme." IEEE Journal of Solid-State Circuits 56.1 (2020): 199-211.
- **[2]** Reed, Irving S., and Gustave Solomon. "Polynomial codes over certain finite fields." Journal of the society for industrial and applied mathematics 8.2 (1960): 300-304.
