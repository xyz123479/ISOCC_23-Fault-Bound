/*
** 체크해야 할 것
 (1) conservative : BL32 묶음 내에서 error가 다른 chip에서 나오는가 (1-bit error는 내부 OECC에서 고쳐서 나온다.)
   => on/off 가능하도록 만들기 [일단 끄고 하자.]
 (2) 

** CE/DUE/SDC
 (1) if Syndrome == 0
  -> NE 'or' SDC
  -> error correction 안하고 SDC check 하면 된다.
 (2) else if Syndrome != 0
  [1] CE가 되는 case면 CE

문제: Correction capability 넘어가는 error 발생시 해당 error가 DUE인지 SDC인지 구분하는 것이 문제

*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctime>
#include<iostream>
#include<cstdlib>
#include <vector>
#include <set>
#include <algorithm>
#include <math.h>
#include <cstring>

// Configuration 
// 필요로 하면 변경하기
#define CHANNEL_WIDTH 40
#define CHIP_NUM 10
#define DATA_CHIP_NUM 8 // sub-channel마다 data chip 개수
#define CHIP_WIDTH 4
#define BLHEIGHT 32 // RECC가 검사해야 하는 Burst length (ondie ecc의 redundancy는 제외하고 BL32만 검사한다. conservative mode 조건 고려하기)

#define OECC_CW_LEN 136 // ondie ecc의 codeword 길이 (bit 단위)
#define OECC_DATA_LEN 128 // ondie ecc의 dataward 길이 (bit 단위)
#define OECC_REDUN_LEN 8 // ondie ecc의 redundancy 길이 (bit 단위)

#define RUN_NUM 100000000 // 실행 횟수


int fault_cnt=0;

//configuration over

using namespace std;
int pp_8 [9] = {1, 0, 1, 1, 1, 0, 0, 0, 1}; /* specify irreducible polynomial coeffts. ex: GF(2^8) 1+x^2+x^3+x^4+x^8 = {1, 0, 1, 1, 1, 0, 0, 0, 1} */ 
int pp_16 [17] = {1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}; /* GF(2^16) 1+x^2+x^3+x^5+x^16 = {1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1} */ 
int alpha_to_8 [256], index_of_8 [256];
int alpha_to_16 [65536], index_of_16 [65536];

unsigned int H_Matrix_SEC_16b_bound[OECC_REDUN_LEN][OECC_CW_LEN]; // 8 x 136
unsigned int H_Matrix_SEC_32b_bound[OECC_REDUN_LEN][OECC_CW_LEN]; // 8 x 136
unsigned int H_Matrix_SEC_Un_bound[OECC_REDUN_LEN][OECC_CW_LEN]; // 8 x 136
enum OECC_TYPE {Bound_1_16=0, Bound_2_8=1, Bound_4_4=2, Bound_1_32=3, Bound_2_16=4, Bound_4_8=5, Un_bound_1_16=6, Un_bound_2_8=7, Un_bound_4_4=8, Un_bound_1_32=9, Un_bound_2_16=10, Un_bound_4_8=11}; // oecc_type
enum FAULT_TYPE {SE_SE=0, SE_MBBE=1, SE_SW=2, SE_SP=3, SE_CHIPKILL=4, MBBE_MBBE=5, MBBE_SW=6, MBBE_SP=7, MBBE_CHIPKILL=8, SW_SW=9, SW_SP=10, SW_CHIPKILL=11, SP_SP=12, SP_CHIPKILL=13, CHIPKILL_CHIPKILL=14}; // fault_type (Error)
enum RECC_TYPE {Chipkill_correct_4_2=0, Chipkill_correct_4_4=1, Chipkill_correct_2_4=2, Chipkill_correct_2_8=3, Chipkill_correct_1_8=4, Chipkill_correct_1_16=5, NO_RLECC=6}; // recc_type
enum RESULT_TYPE {NE=0, CE=1, DUE=2, SDC=3}; // result_type

/* generate GF(2**mm) from the irreducible polynomial p(X) in pp[0]..pp[mm]
   lookup tables:  index->polynomial form   alpha_to[] contains j=alpha**i;
                   polynomial form -> index form  index_of[j=alpha**i] = i
   alpha=2 is the primitive element of GF(2**mm)

    GF(2^4) 기준 (x^4+x+1)

    index_of[0] = -1 -> 0000은 a^-1로 처리
    index_of[a^0 (0001)] = 0
    index_of[a^1 (0010)] = 1
    index_of[a^2 (0100)] = 2
    index_of[a^3 (1000)] = 3
    index_of[a^4 (0011)] = 4
    index_of[a^5 (0110)] = 5
    index_of[a^6 (1100)] = 6
    index_of[a^7 (1011)] = 7
    index_of[a^8 (0101)] = 8
    index_of[a^9 (1010)] = 9
    index_of[a^10 (0111)] = 10
    index_of[a^11 (1110)] = 11
    index_of[a^12 (1111)] = 12
    index_of[a^13 (1101)] = 13
    index_of[a^14 (1001)] = 14
*/
void generate_gf()
{
  //register int i, mask;
  int i, mask;
  int mm_8=8;
  int mm_16=16;

  // 1. GF(2^8)
  mask = 1;
  alpha_to_8[mm_8] = 0;
  
  // a^?의 vector form에서 1이 1개인 부분 -> 0000_0001, 0000_0010, ... 1000_0000 GF(2^8) 기준
  for (i=0; i<mm_8; i++){ // 0~3 GF(2^4), 0~7 GF(2^8) 
    alpha_to_8[i] = mask;
    index_of_8[alpha_to_8[i]] = i;
    if (pp_8[i]!=0)
      alpha_to_8[mm_8] ^= mask;
    mask <<= 1;
  }
  index_of_8[alpha_to_8[mm_8]] = mm_8;
  // a^?의 vector form에서 1이 2개 이상인 부분
  mask >>= 1 ;
  for (i=mm_8+1; i<(int)pow(2.0,mm_8)-1; i++){ // 5~14, 9~254  (원래 nn)
    if (alpha_to_8[i-1] >= mask)
      alpha_to_8[i] = alpha_to_8[mm_8] ^ ((alpha_to_8[i-1]^mask)<<1);
    else 
      alpha_to_8[i] = alpha_to_8[i-1]<<1;
    index_of_8[alpha_to_8[i]] = i;
  }
  index_of_8[0] = -1;

  // 2. GF(2^16)
  mask = 1;
  alpha_to_16[mm_16] = 0;
  
  for (i=0; i<mm_16; i++){ 
    alpha_to_16[i] = mask;
    index_of_16[alpha_to_16[i]] = i;
    if (pp_16[i]!=0)
      alpha_to_16[mm_16] ^= mask;
    mask <<= 1;
  }
  index_of_16[alpha_to_16[mm_16]] = mm_16;

  mask >>= 1 ;
  for (i=mm_16+1; i<(int)pow(2.0,mm_16)-1; i++){
    if (alpha_to_16[i-1] >= mask)
      alpha_to_16[i] = alpha_to_16[mm_16] ^ ((alpha_to_16[i-1]^mask)<<1);
    else 
      alpha_to_16[i] = alpha_to_16[i-1]<<1;
    index_of_16[alpha_to_16[i]] = i;
  }
  index_of_16[0] = -1;
}

// O-ECC H-matrix 생성
void generator_oecc_H_matrix()
{
    FILE *fp5=fopen("H_Matrix_SEC_16b_bound.txt","r");
    while(1){
        unsigned int value;
        for(int row=0; row<OECC_REDUN_LEN; row++){
            for(int column=0; column<OECC_CW_LEN; column++){
                fscanf(fp5,"%d ",&value);
                H_Matrix_SEC_16b_bound[row][column]=value;
                //printf("%d ",H_Matrix_binary[row][column]);
            }
        }
        if(feof(fp5))
            break;
    }
    fclose(fp5);

    FILE *fp6=fopen("H_Matrix_SEC_32b_bound.txt","r");
    while(1){
        unsigned int value;
        for(int row=0; row<OECC_REDUN_LEN; row++){
            for(int column=0; column<OECC_CW_LEN; column++){
                fscanf(fp6,"%d ",&value);
                H_Matrix_SEC_32b_bound[row][column]=value;
                //printf("%d ",H_Matrix_binary[row][column]);
            }
        }
        if(feof(fp6))
            break;
    }
    fclose(fp6);

    FILE *fp4=fopen("H_Matrix_SEC_Un_bound.txt","r");
    while(1){
        unsigned int value;
        for(int row=0; row<OECC_REDUN_LEN; row++){
            for(int column=0; column<OECC_CW_LEN; column++){
                fscanf(fp4,"%d ",&value);
                H_Matrix_SEC_Un_bound[row][column]=value;
                //printf("%d ",H_Matrix_binary[row][column]);
            }
        }
        if(feof(fp4))
            break;
    }
    fclose(fp4);

    return;
}

// OECC, RECC, FAULT TYPE 각각의 type을 string으로 지정. 이것을 기준으로 뒤에서 error injection, oecc, recc에서 어떻게 할지 바뀐다!!!
void oecc_recc_fault_type_assignment(string &OECC, string &FAULT, string &RECC, int *oecc_type, int *fault_type, int*recc_type, int oecc, int fault, int recc)
{
    // 1. OECC TYPE 지정
    // int oecc, int fault, int recc는 main함수 매개변수 argv로 받은 것이다. run.py에서 변경 가능
    switch (oecc){
        case Bound_1_16:
            OECC.replace(OECC.begin(), OECC.end(),"Bound_1_16");
            *oecc_type=Bound_1_16;
            break;
        case Bound_2_8:
            OECC.replace(OECC.begin(), OECC.end(),"Bound_2_8");
            *oecc_type=Bound_2_8;
            break;
        case Bound_4_4:
            OECC.replace(OECC.begin(), OECC.end(),"Bound_4_4");
            *oecc_type=Bound_4_4;
            break;
        case Bound_1_32:
            OECC.replace(OECC.begin(), OECC.end(),"Bound_1_32");
            *oecc_type=Bound_1_32;
            break;
        case Bound_2_16:
            OECC.replace(OECC.begin(), OECC.end(),"Bound_2_16");
            *oecc_type=Bound_2_16;
            break;
        case Bound_4_8:
            OECC.replace(OECC.begin(), OECC.end(),"Bound_4_8");
            *oecc_type=Bound_4_8;
            break;
        case Un_bound_1_16:
            OECC.replace(OECC.begin(), OECC.end(),"Un_bound_1_16");
            *oecc_type=Un_bound_1_16;
            break;
        case Un_bound_2_8:
            OECC.replace(OECC.begin(), OECC.end(),"Un_bound_2_8");
            *oecc_type=Un_bound_2_8;
            break;
        case Un_bound_4_4:
            OECC.replace(OECC.begin(), OECC.end(),"Un_bound_4_4");
            *oecc_type=Un_bound_4_4;
            break;
        case Un_bound_1_32:
            OECC.replace(OECC.begin(), OECC.end(),"Un_bound_1_32");
            *oecc_type=Un_bound_1_32;
            break;
        case Un_bound_2_16:
            OECC.replace(OECC.begin(), OECC.end(),"Un_bound_2_16");
            *oecc_type=Un_bound_2_16;
            break;
        case Un_bound_4_8:
            OECC.replace(OECC.begin(), OECC.end(),"Un_bound_4_8");
            *oecc_type=Un_bound_4_8;
            break;
        default:
            break;
    }
    switch (fault){
        case SE_SE:
            FAULT.replace(FAULT.begin(), FAULT.end(),"SE_SE");
            *fault_type=SE_SE;
            break;
        case SE_MBBE:
            FAULT.replace(FAULT.begin(), FAULT.end(),"SE_MBBE");
            *fault_type=SE_MBBE;
            break;
        case SE_SW:
            FAULT.replace(FAULT.begin(), FAULT.end(),"SE_SW");
            *fault_type=SE_SW;
            break;
        case SE_SP:
            FAULT.replace(FAULT.begin(), FAULT.end(),"SE_SP");
            *fault_type=SE_SP;
            break;
        case SE_CHIPKILL:
            FAULT.replace(FAULT.begin(), FAULT.end(),"SE_CHIPKILL");
            *fault_type=SE_CHIPKILL;
            break;
        case MBBE_MBBE:
            FAULT.replace(FAULT.begin(), FAULT.end(),"MBBE_MBBE");
            *fault_type=MBBE_MBBE;
            break;
        case MBBE_SW:
            FAULT.replace(FAULT.begin(), FAULT.end(),"MBBE_SW");
            *fault_type=MBBE_SW;
            break;
        case MBBE_SP:
            FAULT.replace(FAULT.begin(), FAULT.end(),"MBBE_SP");
            *fault_type=MBBE_SP;
            break;
        case MBBE_CHIPKILL:
            FAULT.replace(FAULT.begin(), FAULT.end(),"MBBE_CHIPKILL");
            *fault_type=MBBE_CHIPKILL;
            break;
        case SW_SW:
            FAULT.replace(FAULT.begin(), FAULT.end(),"SW_SW");
            *fault_type=SW_SW;
            break;
        case SW_SP:
            FAULT.replace(FAULT.begin(), FAULT.end(),"SW_SP");
            *fault_type=SW_SP;
            break;
        case SW_CHIPKILL:
            FAULT.replace(FAULT.begin(), FAULT.end(),"SW_CHIPKILL");
            *fault_type=SW_CHIPKILL;
            break;
        case SP_SP:
            FAULT.replace(FAULT.begin(), FAULT.end(),"SP_SP");
            *fault_type=SP_SP;
            break;
        case SP_CHIPKILL:
            FAULT.replace(FAULT.begin(), FAULT.end(),"SP_CHIPKILL");
            *fault_type=SP_CHIPKILL;
            break;
        case CHIPKILL_CHIPKILL:
            FAULT.replace(FAULT.begin(), FAULT.end(),"CHIPKILL_CHIPKILL");
            *fault_type=CHIPKILL_CHIPKILL;
            break;
        default:
            break;
    }
    switch (recc){
        case Chipkill_correct_4_2: 
            RECC.replace(RECC.begin(), RECC.end(),"Chipkill_correct_4_2");
            *recc_type=Chipkill_correct_4_2;
            break;
        case Chipkill_correct_4_4: 
            RECC.replace(RECC.begin(), RECC.end(),"Chipkill_correct_4_4");
            *recc_type=Chipkill_correct_4_4;
            break;
        case Chipkill_correct_2_4: 
            RECC.replace(RECC.begin(), RECC.end(),"Chipkill_correct_2_4");
            *recc_type=Chipkill_correct_2_4;
            break;
        case Chipkill_correct_2_8: 
            RECC.replace(RECC.begin(), RECC.end(),"Chipkill_correct_2_8");
            *recc_type=Chipkill_correct_2_8;
            break;
        case Chipkill_correct_1_8: 
            RECC.replace(RECC.begin(), RECC.end(),"Chipkill_correct_1_8");
            *recc_type=Chipkill_correct_1_8;
            break;
        case Chipkill_correct_1_16: 
            RECC.replace(RECC.begin(), RECC.end(),"Chipkill_correct_1_16");
            *recc_type=Chipkill_correct_1_16;
            break;
        case NO_RLECC:
            RECC.replace(RECC.begin(), RECC.end(),"NO_RLECC");
            *recc_type=NO_RLECC;
            break;
        default:
            break;
    }
    return;
}

// SE injection (Single Error injection)
void error_injection_SE(int Fault_Chip_position, unsigned int Chip_array[][OECC_CW_LEN], int oecc_type)
{
    int Fault_pos = rand()%OECC_CW_LEN; // 0~135;
    Chip_array[Fault_Chip_position][Fault_pos]^=1;
    return;
}

// MBE injection (Multi bit Bounded Error injection) (More than 1 bit)
void error_injection_MBBE(int Fault_Chip_position, unsigned int Chip_array[][OECC_CW_LEN], int oecc_type)
{
    int error_bound_pos;
    switch(oecc_type){
        case Bound_1_16:
        case Un_bound_1_16:{
            error_bound_pos=rand()%8; // 0~7
            while(1){
                int error_count=0;
                for(int index=0; index<16; index++){
                    Chip_array[Fault_Chip_position][(error_bound_pos/4)*64+index*4+error_bound_pos%4]=rand()%2; // 0 (no-error) 'or' 1 (error)
                    error_count+=Chip_array[Fault_Chip_position][(error_bound_pos/4)*64+index*4+error_bound_pos%4];
                }
                // # of errors: more than 1-bit
                if(error_count>=2)
                    break;
            }
            }
            break;
        case Bound_2_8:
        case Un_bound_2_8:{
            error_bound_pos=rand()%8; // 0~7
            while(1){
                int error_count=0;
                for(int index=0; index<8; index++){
                    Chip_array[Fault_Chip_position][(error_bound_pos/2)*32+index*4+(error_bound_pos%2)*2]=rand()%2; // 0 (no-error) 'or' 1 (error)
                    Chip_array[Fault_Chip_position][(error_bound_pos/2)*32+index*4+(error_bound_pos%2)*2+1]=rand()%2; // 0 (no-error) 'or' 1 (error)
                    error_count+=Chip_array[Fault_Chip_position][(error_bound_pos/2)*32+index*4+(error_bound_pos%2)*2];
                    error_count+=Chip_array[Fault_Chip_position][(error_bound_pos/2)*32+index*4+(error_bound_pos%2)*2+1];
                }
                // # of errors: more than 1-bit
                if(error_count>=2)
                    break;
            }
            }
            break;
        case Bound_4_4:
        case Un_bound_4_4:{
            error_bound_pos=rand()%8; // 0~7
            while(1){
                int error_count=0;
                for(int index=0; index<4; index++){
                    Chip_array[Fault_Chip_position][error_bound_pos*16+index*4]=rand()%2; // 0 (no-error) 'or' 1 (error)
                    Chip_array[Fault_Chip_position][error_bound_pos*16+index*4+1]=rand()%2; // 0 (no-error) 'or' 1 (error)
                    Chip_array[Fault_Chip_position][error_bound_pos*16+index*4+2]=rand()%2; // 0 (no-error) 'or' 1 (error)
                    Chip_array[Fault_Chip_position][error_bound_pos*16+index*4+3]=rand()%2; // 0 (no-error) 'or' 1 (error)
                    error_count+=Chip_array[Fault_Chip_position][error_bound_pos*16+index*4];
                    error_count+=Chip_array[Fault_Chip_position][error_bound_pos*16+index*4+1];
                    error_count+=Chip_array[Fault_Chip_position][error_bound_pos*16+index*4+2];
                    error_count+=Chip_array[Fault_Chip_position][error_bound_pos*16+index*4+3];
                }
                // # of errors: more than 1-bit
                if(error_count>=2)
                    break;
            }
            }
            break;
        case Bound_1_32:
        case Un_bound_1_32:{
            error_bound_pos=rand()%4; // 0~3
            while(1){
                int error_count=0;
                for(int index=0; index<32; index++){
                    Chip_array[Fault_Chip_position][index*4+error_bound_pos%4]=rand()%2; // 0 (no-error) 'or' 1 (error)
                    error_count+=Chip_array[Fault_Chip_position][index*4+error_bound_pos%4];
                }
                // # of errors: more than 1-bit
                if(error_count>=2)
                    break;
            }
            }
            break;
        case Bound_2_16:
        case Un_bound_2_16:{
            error_bound_pos=rand()%4; // 0~3
            while(1){
                int error_count=0;
                for(int index=0; index<16; index++){
                    Chip_array[Fault_Chip_position][(error_bound_pos/2)*64+index*4+(error_bound_pos%2)*2]=rand()%2; // 0 (no-error) 'or' 1 (error)
                    Chip_array[Fault_Chip_position][(error_bound_pos/2)*64+index*4+(error_bound_pos%2)*2+1]=rand()%2; // 0 (no-error) 'or' 1 (error)
                    error_count+=Chip_array[Fault_Chip_position][(error_bound_pos/2)*64+index*4+(error_bound_pos%2)*2];
                    error_count+=Chip_array[Fault_Chip_position][(error_bound_pos/2)*64+index*4+(error_bound_pos%2)*2+1];
                }
                // # of errors: more than 1-bit
                if(error_count>=2)
                    break;
            }
            }
            break;
        case Bound_4_8:
        case Un_bound_4_8:{
            error_bound_pos=rand()%4; // 0~3
            while(1){
                int error_count=0;
                for(int index=0; index<4; index++){
                    Chip_array[Fault_Chip_position][error_bound_pos*32+index*4]=rand()%2; // 0 (no-error) 'or' 1 (error)
                    Chip_array[Fault_Chip_position][error_bound_pos*32+index*4+1]=rand()%2; // 0 (no-error) 'or' 1 (error)
                    Chip_array[Fault_Chip_position][error_bound_pos*32+index*4+2]=rand()%2; // 0 (no-error) 'or' 1 (error)
                    Chip_array[Fault_Chip_position][error_bound_pos*32+index*4+3]=rand()%2; // 0 (no-error) 'or' 1 (error)
                    error_count+=Chip_array[Fault_Chip_position][error_bound_pos*32+index*4];
                    error_count+=Chip_array[Fault_Chip_position][error_bound_pos*32+index*4+1];
                    error_count+=Chip_array[Fault_Chip_position][error_bound_pos*32+index*4+2];
                    error_count+=Chip_array[Fault_Chip_position][error_bound_pos*32+index*4+3];
                }
                // # of errors: more than 1-bit
                if(error_count>=2)
                    break;
            }
            }
            break;
        default:
            break;
    }
    return;
}

// SW injection (Single Word Error injection) -> All 4 bits are corrupted
void error_injection_SW(int Fault_Chip_position, unsigned int Chip_array[][OECC_CW_LEN], int oecc_type)
{
    int fault_beat_pos = rand()%(OECC_CW_LEN/CHIP_WIDTH); // 0~33
    for(int index=0; index<CHIP_WIDTH; index++) // 0~3
        Chip_array[Fault_Chip_position][fault_beat_pos*CHIP_WIDTH+index]^=1;

    return;
}

// SP injection (Single Pin Error injection) -> more than 2 bits
// Parity에서는 error 발생 X
void error_injection_SP(int Fault_Chip_position, unsigned int Chip_array[][OECC_CW_LEN], int oecc_type)
{
    while(1){
        // 0으로 초기화
        for(int index=0; index<OECC_DATA_LEN; index++)
            Chip_array[Fault_Chip_position][index]=0;
        
        // error injection
        int Fault_pin_pos=rand()%CHIP_WIDTH; // 0~3
        int pin_error_pos=Fault_pin_pos;
        while(pin_error_pos<OECC_DATA_LEN){
            if(rand()%2!=0) // 0(no error) 'or' 1(bit-flip)
                Chip_array[Fault_Chip_position][pin_error_pos]^=1;
            pin_error_pos+=CHIP_WIDTH;
        }

        // 2bit 초과로 error가 발생했는지 check
        int count=0;
        for(int index=0; index<OECC_DATA_LEN; index++){
            if(Chip_array[Fault_Chip_position][index]==1){count++;}
        }
        if(count>=3)
            break;
    }
    return;
}

// Chipkill injection
void error_injection_CHIPKILL(int Fault_Chip_position, unsigned int Chip_array[][OECC_CW_LEN], int oecc_type)
{
    for(int Fault_pos=0; Fault_pos<OECC_CW_LEN; Fault_pos++){ // 0~135
        if(rand()%2!=0) // 0(no error) 'or' 1(bit-flip)
            Chip_array[Fault_Chip_position][Fault_pos]^=1;
    }

    return;
}

// OECC 1bit correction
void error_correction_oecc(int Fault_Chip_position, unsigned int Chip_array[][OECC_CW_LEN], int oecc_type)
{
    unsigned int Syndromes[OECC_REDUN_LEN]; // 8 x 1
    int cnt=0;

    switch (oecc_type){
        case Bound_1_16:
            // Syndromes = H * C^T
            for(int row=0; row<OECC_REDUN_LEN; row++){
                unsigned int row_value=0;
                for(int column=0; column<16; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][60-column*CHIP_WIDTH]); // 60, 56, 52, ... 0 -> 총 16개
                for(int column=16; column<32; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][61-(column-16)*CHIP_WIDTH]); // 61, 57, 53, ... 1 -> 총 16개
                for(int column=32; column<48; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][62-(column-32)*CHIP_WIDTH]); // 62, 58, 54, ... 2 -> 총 16개
                for(int column=48; column<64; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][63-(column-48)*CHIP_WIDTH]); // 63, 59, 55, ... 3 -> 총 16개
                for(int column=64; column<80; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][124-(column-64)*CHIP_WIDTH]); // 124, 120, 116, ... 64 -> 총 16개
                for(int column=80; column<96; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][125-(column-80)*CHIP_WIDTH]); // 125, 121, 117, ... 65 -> 총 16개
                for(int column=96; column<112; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][126-(column-96)*CHIP_WIDTH]); // 126, 122, 118, ... 66 -> 총 16개
                for(int column=112; column<128; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][127-(column-112)*CHIP_WIDTH]); // 127, 123, 119, ... 67 -> 총 16개
                for(int column=128; column<136; column++) // parity 부분
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][column]); // 128~135 -> 총 8개
                Syndromes[row]=row_value;
            }
            // error correction (Check Syndromes)
            for(int error_pos=0; error_pos<OECC_CW_LEN; error_pos++){
                cnt=0;
                for(int row=0; row<OECC_REDUN_LEN; row++){
                    if(Syndromes[row]==H_Matrix_SEC_16b_bound[row][error_pos])
                        cnt++;
                    else
                        break;
                }
                // 1-bit error 일때만 error correction 진행
                if(cnt==OECC_REDUN_LEN){
                    if(error_pos<16) // error_pos:0 => 60, error_pos:1 => 56, ... error_pos:15 => 0
                        Chip_array[Fault_Chip_position][60-error_pos*CHIP_WIDTH]^=1; // error correction (bit-flip)
                    else if(error_pos<32) // error_pos:16 => 61, error_pos:17 => 57, ... error_pos:31 => 1
                        Chip_array[Fault_Chip_position][61-(error_pos-16)*CHIP_WIDTH]^=1; // error correction (bit-flip)
                    else if(error_pos<48) 
                        Chip_array[Fault_Chip_position][62-(error_pos-32)*CHIP_WIDTH]^=1; // error correction (bit-flip)
                    else if(error_pos<64) 
                        Chip_array[Fault_Chip_position][63-(error_pos-48)*CHIP_WIDTH]^=1; // error correction (bit-flip)
                    else if(error_pos<80) //
                        Chip_array[Fault_Chip_position][124-(error_pos-64)*CHIP_WIDTH]^=1; // error correction (bit-flip)
                    else if(error_pos<96) 
                        Chip_array[Fault_Chip_position][125-(error_pos-80)*CHIP_WIDTH]^=1; // error correction (bit-flip)
                    else if(error_pos<112) 
                        Chip_array[Fault_Chip_position][126-(error_pos-96)*CHIP_WIDTH]^=1; // error correction (bit-flip)
                    else if(error_pos<128) 
                        Chip_array[Fault_Chip_position][127-(error_pos-112)*CHIP_WIDTH]^=1; // error correction (bit-flip)
                    else if(error_pos<136)
                        Chip_array[Fault_Chip_position][error_pos]^=1; // error correction (bit-flip)
                }
            }
            break;
        case Bound_2_8:
            // Syndromes = H * C^T
            for(int row=0; row<OECC_REDUN_LEN; row++){
                unsigned int row_value=0;
                for(int column=0; column<16; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][60-4*(column/2)+(column%2)]); // 60, 61, 56, 57, 52, 53, 48, 49, 44, 45, 40, 41, 36, 37, 32, 33 -> 총 16개
                for(int column=16; column<32; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][62 - 4*((column-16)/2)+((column-16)%2)]); // 62, 63, 58, 59, 54, 55, 50, 51, 46, 47, 42, 43, 38, 39, 34, 35 -> 총 16개
                for(int column=32; column<48; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][28 - 4*((column-32)/2)+((column-32)%2)]); // 28, 29, 24, 25, 20, 21, 16, 17, 12, 13, 8, 9, 4, 5, 0, 1 -> 총 16개
                for(int column=48; column<64; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][30 - 4*((column-48)/2)+((column-48)%2)]); // 30, 31, 26, 27, 22, 23, 18, 19, 14, 15, 10, 11, 6, 7, 2, 3 -> 총 16개
                for(int column=64; column<80; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][124 - 4*((column-64)/2)+((column-64)%2)]); // 124, 125, 120, 121, 116, 117, 112, 113, 108, 109, 104, 105, 100, 101, 96, 97 -> 총 16개
                for(int column=80; column<96; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][126 - 4*((column-80)/2)+((column-80)%2)]); // 126, 127, 122, 123, 118, 119, 114, 115, 110, 111, 106, 107, 102, 103, 98, 99 -> 총 16개
                for(int column=96; column<112; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][92 - 4*((column-96)/2)+((column-96)%2)]); // 92, 93, 88, 89, 84, 85, 80, 81, 76, 77, 72, 73, 68, 69, 64, 65 -> 총 16개
                for(int column=112; column<128; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][94 - 4*((column-112)/2)+((column-112)%2)]); // 94, 95, 90, 91, 86, 87, 82, 83, 78, 79, 74, 75, 70, 71, 66, 67 -> 총 16개
                for(int column=128; column<136; column++) // parity 부분
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][column]); // 128~135 -> 총 8개
                Syndromes[row]=row_value;
            }
            // error correction (Check Syndromes)
            for(int error_pos=0; error_pos<OECC_CW_LEN; error_pos++){
                cnt=0;
                for(int row=0; row<OECC_REDUN_LEN; row++){
                    if(Syndromes[row]==H_Matrix_SEC_16b_bound[row][error_pos])
                        cnt++;
                    else
                        break;
                }
                // 1-bit error 일때만 error correction 진행
                if(cnt==OECC_REDUN_LEN){
                    if(error_pos<16) 
                        Chip_array[Fault_Chip_position][60-4*(error_pos/2)+(error_pos%2)]^=1; // error correction (bit-flip)
                    else if(error_pos<32) 
                        Chip_array[Fault_Chip_position][62 - 4*((error_pos-16)/2)+((error_pos-16)%2)]^=1; // error correction (bit-flip)
                    else if(error_pos<48) 
                        Chip_array[Fault_Chip_position][28 - 4*((error_pos-32)/2)+((error_pos-32)%2)]^=1; // error correction (bit-flip)
                    else if(error_pos<64) 
                        Chip_array[Fault_Chip_position][30 - 4*((error_pos-48)/2)+((error_pos-48)%2)]^=1; // error correction (bit-flip)
                    else if(error_pos<80) //
                        Chip_array[Fault_Chip_position][124 - 4*((error_pos-64)/2)+((error_pos-64)%2)]^=1; // error correction (bit-flip)
                    else if(error_pos<96) 
                        Chip_array[Fault_Chip_position][126 - 4*((error_pos-80)/2)+((error_pos-80)%2)]^=1; // error correction (bit-flip)
                    else if(error_pos<112) 
                        Chip_array[Fault_Chip_position][92 - 4*((error_pos-96)/2)+((error_pos-96)%2)]^=1; // error correction (bit-flip)
                    else if(error_pos<128) 
                        Chip_array[Fault_Chip_position][94 - 4*((error_pos-112)/2)+((error_pos-112)%2)]^=1; // error correction (bit-flip)
                    else if(error_pos<136)
                        Chip_array[Fault_Chip_position][error_pos]^=1; // error correction (bit-flip)
                }
            }
            break;
        case Bound_4_4:
            // Syndromes = H * C^T
            for(int row=0; row<OECC_REDUN_LEN; row++){
                unsigned int row_value=0;
                for(int column=0; column<16; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][60-4*(column/4)+(column%4)]); // 60, 61, 62, 63, 56, 57, 58, 59, 52, 53, 54, 55, 48, 49, 50, 51 -> 총 16개
                for(int column=16; column<32; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][44-4*((column-16)/4)+((column-16)%4)]); // 44, 45, 46, 47, 40, 41, 42, 43, 36, 37, 38, 39, 32, 33, 34, 35 -> 총 16개
                for(int column=32; column<48; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][28-4*((column-32)/4)+((column-32)%4)]); // 28, 29, 30, 31, 24, 25, 26, 27, 20, 21, 22, 23, 16, 17, 18, 19 -> 총 16개
                for(int column=48; column<64; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][12-4*((column-48)/4)+((column-48)%4)]); // 12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3 -> 총 16개
                for(int column=64; column<80; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][124-4*((column-64)/4)+((column-64)%4)]); 
                for(int column=80; column<96; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][108-4*((column-80)/4)+((column-80)%4)]); 
                for(int column=96; column<112; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][92-4*((column-96)/4)+((column-96)%4)]); 
                for(int column=112; column<128; column++)
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][76-4*((column-112)/4)+((column-112)%4)]); 
                for(int column=128; column<136; column++) // parity 부분
                    row_value=row_value^(H_Matrix_SEC_16b_bound[row][column] * Chip_array[Fault_Chip_position][column]); // 128~135 -> 총 8개
                Syndromes[row]=row_value;
            }
            // error correction (Check Syndromes)
            for(int error_pos=0; error_pos<OECC_CW_LEN; error_pos++){
                cnt=0;
                for(int row=0; row<OECC_REDUN_LEN; row++){
                    if(Syndromes[row]==H_Matrix_SEC_16b_bound[row][error_pos])
                        cnt++;
                    else
                        break;
                }
                // 1-bit error 일때만 error correction 진행
                if(cnt==OECC_REDUN_LEN){
                    if(error_pos<16)
                        Chip_array[Fault_Chip_position][60-4*(error_pos/4)+(error_pos%4)]^=1; // error correction (bit-flip)
                    else if(error_pos<32)
                        Chip_array[Fault_Chip_position][44-4*((error_pos-16)/4)+((error_pos-16)%4)]^=1; // error correction (bit-flip)
                    else if(error_pos<48)
                        Chip_array[Fault_Chip_position][28-4*((error_pos-32)/4)+((error_pos-32)%4)]^=1; // error correction (bit-flip)
                    else if(error_pos<64)
                        Chip_array[Fault_Chip_position][12-4*((error_pos-48)/4)+((error_pos-48)%4)]^=1; // error correction (bit-flip)
                    else if(error_pos<80)
                        Chip_array[Fault_Chip_position][124-4*((error_pos-64)/4)+((error_pos-64)%4)]^=1; // error correction (bit-flip)
                    else if(error_pos<96)
                        Chip_array[Fault_Chip_position][108-4*((error_pos-80)/4)+((error_pos-80)%4)]^=1; // error correction (bit-flip)
                    else if(error_pos<112) 
                        Chip_array[Fault_Chip_position][92-4*((error_pos-96)/4)+((error_pos-96)%4)]^=1; // error correction (bit-flip)
                    else if(error_pos<128)
                        Chip_array[Fault_Chip_position][76-4*((error_pos-112)/4)+((error_pos-112)%4)]^=1; // error correction (bit-flip)
                    else if(error_pos<136)
                        Chip_array[Fault_Chip_position][error_pos]^=1; // error correction (bit-flip)
                }
            }
            break;
        case Bound_1_32:
            // Syndromes = H * C^T
            for(int row=0; row<OECC_REDUN_LEN; row++){
                unsigned int row_value=0;
                for(int column=0; column<32; column++)
                    row_value=row_value^(H_Matrix_SEC_32b_bound[row][column] * Chip_array[Fault_Chip_position][124-4*(column)]); // 124, 120, 116, 112, 108, 104, 100, 96, 92, 88, 84, 80, 76, 72, 68, 64, ..., 0 -> 총 32개
                for(int column=32; column<64; column++)
                    row_value=row_value^(H_Matrix_SEC_32b_bound[row][column] * Chip_array[Fault_Chip_position][125-4*(column-32)]); // 125, 121, 117, 113, ... , 1 -> 총 32개
                for(int column=64; column<96; column++)
                    row_value=row_value^(H_Matrix_SEC_32b_bound[row][column] * Chip_array[Fault_Chip_position][126-4*(column-64)]); // 126, 122, ..., 2 -> 총 32개
                for(int column=96; column<128; column++)
                    row_value=row_value^(H_Matrix_SEC_32b_bound[row][column] * Chip_array[Fault_Chip_position][127-4*(column-96)]); 
                for(int column=128; column<136; column++) // parity 부분
                    row_value=row_value^(H_Matrix_SEC_32b_bound[row][column] * Chip_array[Fault_Chip_position][column]); // 128~135 -> 총 8개
                Syndromes[row]=row_value;
            }
            // error correction (Check Syndromes)
            for(int error_pos=0; error_pos<OECC_CW_LEN; error_pos++){
                cnt=0;
                for(int row=0; row<OECC_REDUN_LEN; row++){
                    if(Syndromes[row]==H_Matrix_SEC_32b_bound[row][error_pos])
                        cnt++;
                    else
                        break;
                }
                // 1-bit error 일때만 error correction 진행
                if(cnt==OECC_REDUN_LEN){
                    if(error_pos<32) // error_pos : 0 => 124, 1 => 120, 2 => 116, ...
                        Chip_array[Fault_Chip_position][124-4*(error_pos)]^=1; // error correction (bit-flip)
                    else if(error_pos<64) // error_pos : 32 => 125, 33 => 121, 34 => 117, .. 
                        Chip_array[Fault_Chip_position][125-4*(error_pos-32)]^=1; // error correction (bit-flip)
                    else if(error_pos<96) // 
                        Chip_array[Fault_Chip_position][126-4*(error_pos-64)]^=1; // error correction (bit-flip)
                    else if(error_pos<128)
                        Chip_array[Fault_Chip_position][127-4*(error_pos-96)]^=1; // error correction (bit-flip)
                    else if(error_pos<136)
                        Chip_array[Fault_Chip_position][error_pos]^=1; // error correction (bit-flip)
                }
            }
            break;
        case Bound_2_16:
            // Syndromes = H * C^T
            for(int row=0; row<OECC_REDUN_LEN; row++){
                unsigned int row_value=0;
                for(int column=0; column<32; column++)
                    row_value=row_value^(H_Matrix_SEC_32b_bound[row][column] * Chip_array[Fault_Chip_position][60-4*(column/2)+(column%2)]); // 60, 61, 56, 57, ... 0, 1 -> 총 32개
                for(int column=32; column<64; column++)
                    row_value=row_value^(H_Matrix_SEC_32b_bound[row][column] * Chip_array[Fault_Chip_position][62-4*((column-32)/2)+((column-32)%2)]); 
                for(int column=64; column<96; column++)
                    row_value=row_value^(H_Matrix_SEC_32b_bound[row][column] * Chip_array[Fault_Chip_position][124-4*((column-64)/2)+((column-64)%2)]); 
                for(int column=96; column<128; column++)
                    row_value=row_value^(H_Matrix_SEC_32b_bound[row][column] * Chip_array[Fault_Chip_position][126-4*((column-96)/2)+((column-96)%2)]); 
                for(int column=128; column<136; column++) // parity 부분
                    row_value=row_value^(H_Matrix_SEC_32b_bound[row][column] * Chip_array[Fault_Chip_position][column]); // 128~135 -> 총 8개
                Syndromes[row]=row_value;
            }
            // error correction (Check Syndromes)
            for(int error_pos=0; error_pos<OECC_CW_LEN; error_pos++){
                cnt=0;
                for(int row=0; row<OECC_REDUN_LEN; row++){
                    if(Syndromes[row]==H_Matrix_SEC_32b_bound[row][error_pos])
                        cnt++;
                    else
                        break;
                }
                // 1-bit error 일때만 error correction 진행
                if(cnt==OECC_REDUN_LEN){
                    if(error_pos<32)
                        Chip_array[Fault_Chip_position][60-4*(error_pos/2)+(error_pos%2)]^=1; // error correction (bit-flip)
                    else if(error_pos<64)
                        Chip_array[Fault_Chip_position][62-4*((error_pos-32)/2)+((error_pos-32)%2)]^=1; // error correction (bit-flip)
                    else if(error_pos<96)
                        Chip_array[Fault_Chip_position][124-4*((error_pos-64)/2)+((error_pos-64)%2)]^=1; // error correction (bit-flip)
                    else if(error_pos<128)
                        Chip_array[Fault_Chip_position][126-4*((error_pos-96)/2)+((error_pos-96)%2)]^=1; // error correction (bit-flip)
                    else if(error_pos<136)
                        Chip_array[Fault_Chip_position][error_pos]^=1; // error correction (bit-flip)
                }
            }
            break;
        case Bound_4_8:
            // Syndromes = H * C^T
            for(int row=0; row<OECC_REDUN_LEN; row++){
                unsigned int row_value=0;
                for(int column=0; column<32; column++)
                    row_value=row_value^(H_Matrix_SEC_32b_bound[row][column] * Chip_array[Fault_Chip_position][60-4*(column/4)+(column%4)]); // 60, 61, 62, 63, 56, 57, ..., 32, 33, 34, 35  -> 총 32개
                for(int column=32; column<64; column++)
                    row_value=row_value^(H_Matrix_SEC_32b_bound[row][column] * Chip_array[Fault_Chip_position][28-4*((column-32)/4)+((column-32)%4)]); // 28, 29, 30, 31, 24, 25, ..., 0, 1, 2, 3  -> 총 32개
                for(int column=64; column<96; column++)
                    row_value=row_value^(H_Matrix_SEC_32b_bound[row][column] * Chip_array[Fault_Chip_position][124-4*((column-64)/4)+((column-64)%4)]); // 124, ... 
                for(int column=96; column<128; column++)
                    row_value=row_value^(H_Matrix_SEC_32b_bound[row][column] * Chip_array[Fault_Chip_position][92-4*((column-96)/4)+((column-96)%4)]); // 92, ...
                for(int column=128; column<136; column++) // parity 부분
                    row_value=row_value^(H_Matrix_SEC_32b_bound[row][column] * Chip_array[Fault_Chip_position][column]); // 128~135 -> 총 8개
                Syndromes[row]=row_value;
            }
            // error correction (Check Syndromes)
            for(int error_pos=0; error_pos<OECC_CW_LEN; error_pos++){
                cnt=0;
                for(int row=0; row<OECC_REDUN_LEN; row++){
                    if(Syndromes[row]==H_Matrix_SEC_32b_bound[row][error_pos])
                        cnt++;
                    else
                        break;
                }
                // 1-bit error 일때만 error correction 진행
                if(cnt==OECC_REDUN_LEN){
                    if(error_pos<32)
                        Chip_array[Fault_Chip_position][60-4*(error_pos/4)+(error_pos%4)]^=1; // error correction (bit-flip)
                    else if(error_pos<64)
                        Chip_array[Fault_Chip_position][28-4*((error_pos-32)/4)+((error_pos-32)%4)]^=1; // error correction (bit-flip)
                    else if(error_pos<96)
                        Chip_array[Fault_Chip_position][124-4*((error_pos-64)/4)+((error_pos-64)%4)]^=1; // error correction (bit-flip)
                    else if(error_pos<128)
                        Chip_array[Fault_Chip_position][92-4*((error_pos-96)/4)+((error_pos-96)%4)]^=1; // error correction (bit-flip)
                    else if(error_pos<136)
                        Chip_array[Fault_Chip_position][error_pos]^=1; // error correction (bit-flip)
                }
            }
            break;
        case Un_bound_1_16:
        case Un_bound_2_8:
        case Un_bound_4_4:
        case Un_bound_1_32:
        case Un_bound_2_16:
        case Un_bound_4_8:
            // Syndromes = H * C^T
            for(int row=0; row<OECC_REDUN_LEN; row++){
                unsigned int row_value=0;
                for(int column=0; column<OECC_CW_LEN; column++)
                    row_value=row_value^(H_Matrix_SEC_Un_bound[row][column] * Chip_array[Fault_Chip_position][column]);
                Syndromes[row]=row_value;
            }
            // error correction (Check Syndromes)
            for(int error_pos=0; error_pos<OECC_CW_LEN; error_pos++){
                cnt=0;
                for(int row=0; row<OECC_REDUN_LEN; row++){
                    if(Syndromes[row]==H_Matrix_SEC_Un_bound[row][error_pos])
                        cnt++;
                    else
                        break;
                }
                // 1-bit error 일때만 error correction 진행
                if(cnt==OECC_REDUN_LEN){
                    Chip_array[Fault_Chip_position][error_pos]^=1; // error correction (bit-flip)
                    return;
                }
            }
            break;
        default:
            break;
    }
    

    // Error가 발생하지 않았거나, 발생했지만 1-bit error는 아닌 경우이다.
    // 이 경우에는 correction을 진행하지 않는다.
    return;
}

// GF(2^8) Rank-level ECC decoding
int decode_rs_8(int recd[], int bb[], int nn, int kk, int tt, int nn_short)
{
   register int i,j,u,q ;
   int elp[nn-kk+2][nn-kk], d[nn-kk+2], l[nn-kk+2], u_lu[nn-kk+2], s[nn-kk+1] ;
   int count=0, syn_error=0, root[tt], loc[tt], z[tt+1], err[nn], reg[tt+1] ;

/* first form the syndromes */
  for (i=1; i<=nn-kk; i++){ 
    s[i] = 0;
      
    for (j=0; j<nn; j++)
      if (recd[j]!=-1)
          s[i] ^= alpha_to_8[(recd[j]+i*j)%nn];      /* recd[j] in index form */
/* convert syndrome from polynomial form to index form  */
      if (s[i]!=0)  syn_error=1;        /* set flag if non-zero syndrome => error */
      s[i] = index_of_8[s[i]];
  }

/* compute the error location polynomial via the Berlekamp iterative algorithm,
   following the terminology of Lin and Costello :
   d[u]: the 'mu'th discrepancy, where u='mu'+1 and 'mu' (the Greek letter!) is the step number
   ranging from -1 to 2*tt (see L&C)
   l[u]: the degree of the elp at that step, and u_l[u] is the difference between the
   step number and the degree of the elp.
*/
 // Syndrome != 0 -> CE 'or' DUE 'or' SDC
  if (syn_error)       /* if errors, try and correct */
  {
/* initialise table entries */
      d[0] = 0 ;           /* index form */
      d[1] = s[1] ;        /* index form */
      elp[0][0] = 0 ;      /* index form */
      elp[1][0] = 1 ;      /* polynomial form */
      for (i=1; i<nn-kk; i++)
        { elp[0][i] = -1 ;   /* index form */
          elp[1][i] = 0 ;   /* polynomial form */
        }
      l[0] = 0 ;
      l[1] = 0 ;
      u_lu[0] = -1 ;
      u_lu[1] = 0 ;
      u = 0 ;

      do
      {
        u++ ;
        if (d[u]==-1)
          { l[u+1] = l[u];
            for (i=0; i<=l[u]; i++)
             {  elp[u+1][i] = elp[u][i] ;
                elp[u][i] = index_of_8[elp[u][i]] ;
             }
          }
        else
/* search for words with greatest u_lu[q] for which d[q]!=0 */
          { q = u-1 ;
            while ((d[q]==-1) && (q>0)) q-- ;
/* have found first non-zero d[q]  */
            if (q>0)
             { j=q ;
               do
               { j-- ;
                 if ((d[j]!=-1) && (u_lu[q]<u_lu[j]))
                   q = j ;
               }while (j>0) ;
             }

/* have now found q such that d[u]!=0 and u_lu[q] is maximum */
/* store degree of new elp polynomial */
            if (l[u]>l[q]+u-q)  l[u+1] = l[u] ;
            else  l[u+1] = l[q]+u-q ;

/* form new elp(x) */
            for (i=0; i<nn-kk; i++)    elp[u+1][i] = 0;
            for (i=0; i<=l[q]; i++)
              if (elp[q][i]!=-1)
                elp[u+1][i+u-q] = alpha_to_8[(d[u]+nn-d[q]+elp[q][i])%nn];
            for (i=0; i<=l[u]; i++)
              { elp[u+1][i] ^= elp[u][i];
                elp[u][i] = index_of_8[elp[u][i]];  /*convert old elp value to index*/
              }
          }
        u_lu[u+1] = u-l[u+1];

/* form (u+1)th discrepancy */
        if (u<nn-kk)    /* no discrepancy computed on last iteration */
          {
            if (s[u+1]!=-1)
                   d[u+1] = alpha_to_8[s[u+1]] ;
            else
              d[u+1] = 0 ;
            for (i=1; i<=l[u+1]; i++)
              if ((s[u+1-i]!=-1) && (elp[u+1][i]!=0))
                d[u+1] ^= alpha_to_8[(s[u+1-i]+index_of_8[elp[u+1][i]])%nn] ;
            d[u+1] = index_of_8[d[u+1]] ;    /* put d[u+1] into index form */
          }
      } while ((u<nn-kk) && (l[u+1]<=tt)) ;

      u++ ;

      ////////////////////////////////////////////////////
      // error correction start!!!!!!!!!!!!!!!!!!!!!!!!!
      ////////////////////////////////////////////////////

      // CE 'or' SDC cases
      if (l[u]<=tt)         /* can correct error */
      {
/* put elp into index form */
         for (i=0; i<=l[u]; i++)   elp[u][i] = index_of_8[elp[u][i]] ;

/* find roots of the error location polynomial -> finding error location */
        for (i=1; i<=l[u]; i++)
          reg[i] = elp[u][i] ;
        
        count = 0 ;
        for (i=1; i<=nn; i++){  
          q = 1 ;
          for (j=1; j<=l[u]; j++)
            if (reg[j]!=-1){ 
              reg[j] = (reg[j]+j)%nn;
              q ^= alpha_to_8[reg[j]];
            }
          if (!q) {        /* store root and error location number indices */
            root[count] = i;
            loc[count] = nn-i; // -> error location
            count++ ;
          }
        }
        
        int CE_cases=1;
        //printf("error location check! (shortened RS code!)\n");
        for(int index=0; index<count; index++){
          if(loc[index]>=nn_short) // except for zero-padding part!!
            CE_cases=0;
        }
        

         // CE 'or' SDC cases
         if (count==l[u] && CE_cases==1){    /* no. roots = degree of elp hence <= tt errors*/
          //printf("CE 'or' SDC cases\n");
          //printf("count : %d\n",count);
          //printf("loc [0] : %d\n",loc[0]);
/* form polynomial z(x) */
          for (i=1; i<=l[u]; i++){        /* Z[0] = 1 always - do not need */
            if ((s[i]!=-1) && (elp[u][i]!=-1))
              z[i] = alpha_to_8[s[i]] ^ alpha_to_8[elp[u][i]];
            else if ((s[i]!=-1) && (elp[u][i]==-1))
              z[i] = alpha_to_8[s[i]] ;
            else if ((s[i]==-1) && (elp[u][i]!=-1))
              z[i] = alpha_to_8[elp[u][i]] ;
            else
              z[i] = 0 ;
            for (j=1; j<i; j++){
              if ((s[j]!=-1) && (elp[u][i-j]!=-1))
                z[i] ^= alpha_to_8[(elp[u][i-j] + s[j])%nn];
            }
            z[i] = index_of_8[z[i]];         /* put into index form */
          }

  /* evaluate errors at locations given by error location numbers loc[i] */
          for (i=0; i<nn; i++){
            err[i] = 0;
            if (recd[i]!=-1)        /* convert recd[] to polynomial form */
              recd[i] = alpha_to_8[recd[i]];
            else  
              recd[i] = 0;
          }
          for (i=0; i<l[u]; i++){    /* compute numerator of error term first */
            err[loc[i]] = 1;       /* accounts for z[0] */
            for (j=1; j<=l[u]; j++){
              if (z[j]!=-1)
                err[loc[i]] ^= alpha_to_8[(z[j]+j*root[i])%nn];
            }
            if (err[loc[i]]!=0){
              err[loc[i]] = index_of_8[err[loc[i]]];
                q = 0;     /* form denominator of error term */
                for (j=0; j<l[u]; j++){
                  if (j!=i)
                    q += index_of_8[1^alpha_to_8[(loc[j]+root[i])%nn]];
                }
                q = q % nn;
                err[loc[i]] = alpha_to_8[(err[loc[i]]-q+nn)%nn];
                recd[loc[i]] ^= err[loc[i]];  /*recd[i] must be in polynomial form */
                //printf("loc [i] : %d\n",loc[i]);
                //printf("err[loc[i]] : %d\n",err[loc[i]]);
            }
          }

          return CE;
          /*
          printf("err (error values) : ");
          for(int index=0; index<nn; index++)
            printf("%d ",err[index]);
          printf("\n");
          */
        }
        // DUE cases
        else{    /* no. roots != degree of elp => >tt errors and cannot solve */
            for (i=0; i<nn; i++){
            err[i] = 0;
            if (recd[i]!=-1)        /* convert recd[] to polynomial form */
              recd[i] = alpha_to_8[recd[i]];
            else  
              recd[i] = 0;
          }
          return DUE;
        }
      }
    // DUE cases
     else{         /* elp has degree has degree >tt hence cannot solve */
        for (i=0; i<nn; i++){
            err[i] = 0;
            if (recd[i]!=-1)        /* convert recd[] to polynomial form */
              recd[i] = alpha_to_8[recd[i]];
            else  
              recd[i] = 0;
          }
       return DUE;
     }
  }
  // NE 'or' SDC cases
  else{       /* no non-zero syndromes => no errors: output received codeword */
        for (i=0; i<nn; i++){
            err[i] = 0;
            if (recd[i]!=-1)        /* convert recd[] to polynomial form */
              recd[i] = alpha_to_8[recd[i]];
            else  
              recd[i] = 0;
        }
    //printf("NE 'or' SDC cases!\n");
    return NE;
  }
}

// GF(2^16) Rank-level ECC decoding
int decode_rs_16(int recd[], int bb[], int nn, int kk, int tt, int nn_short)
{
   register int i,j,u,q ;
   int elp[nn-kk+2][nn-kk], d[nn-kk+2], l[nn-kk+2], u_lu[nn-kk+2], s[nn-kk+1] ;
   int count=0, syn_error=0, root[tt], loc[tt], z[tt+1], err[nn], reg[tt+1] ;

/* first form the syndromes */
  for (i=1; i<=nn-kk; i++){ 
    s[i] = 0;
      
    for (j=0; j<nn; j++)
      if (recd[j]!=-1)
          s[i] ^= alpha_to_16[(recd[j]+i*j)%nn];      /* recd[j] in index form */
/* convert syndrome from polynomial form to index form  */
      if (s[i]!=0)  syn_error=1;        /* set flag if non-zero syndrome => error */
      s[i] = index_of_16[s[i]];
  }

/* compute the error location polynomial via the Berlekamp iterative algorithm,
   following the terminology of Lin and Costello :
   d[u]: the 'mu'th discrepancy, where u='mu'+1 and 'mu' (the Greek letter!) is the step number
   ranging from -1 to 2*tt (see L&C)
   l[u]: the degree of the elp at that step, and u_l[u] is the difference between the
   step number and the degree of the elp.
*/
 // Syndrome != 0 -> CE 'or' DUE 'or' SDC
  if (syn_error)       /* if errors, try and correct */
  {
/* initialise table entries */
      d[0] = 0 ;           /* index form */
      d[1] = s[1] ;        /* index form */
      elp[0][0] = 0 ;      /* index form */
      elp[1][0] = 1 ;      /* polynomial form */
      for (i=1; i<nn-kk; i++)
        { elp[0][i] = -1 ;   /* index form */
          elp[1][i] = 0 ;   /* polynomial form */
        }
      l[0] = 0 ;
      l[1] = 0 ;
      u_lu[0] = -1 ;
      u_lu[1] = 0 ;
      u = 0 ;

      do
      {
        u++ ;
        if (d[u]==-1)
          { l[u+1] = l[u];
            for (i=0; i<=l[u]; i++)
             {  elp[u+1][i] = elp[u][i] ;
                elp[u][i] = index_of_16[elp[u][i]] ;
             }
          }
        else
/* search for words with greatest u_lu[q] for which d[q]!=0 */
          { q = u-1 ;
            while ((d[q]==-1) && (q>0)) q-- ;
/* have found first non-zero d[q]  */
            if (q>0)
             { j=q ;
               do
               { j-- ;
                 if ((d[j]!=-1) && (u_lu[q]<u_lu[j]))
                   q = j ;
               }while (j>0) ;
             }

/* have now found q such that d[u]!=0 and u_lu[q] is maximum */
/* store degree of new elp polynomial */
            if (l[u]>l[q]+u-q)  l[u+1] = l[u] ;
            else  l[u+1] = l[q]+u-q ;

/* form new elp(x) */
            for (i=0; i<nn-kk; i++)    elp[u+1][i] = 0 ;
            for (i=0; i<=l[q]; i++)
              if (elp[q][i]!=-1)
                elp[u+1][i+u-q] = alpha_to_16[(d[u]+nn-d[q]+elp[q][i])%nn] ;
            for (i=0; i<=l[u]; i++)
              { elp[u+1][i] ^= elp[u][i] ;
                elp[u][i] = index_of_16[elp[u][i]] ;  /*convert old elp value to index*/
              }
          }
        u_lu[u+1] = u-l[u+1] ;

/* form (u+1)th discrepancy */
        if (u<nn-kk)    /* no discrepancy computed on last iteration */
          {
            if (s[u+1]!=-1)
                   d[u+1] = alpha_to_16[s[u+1]] ;
            else
              d[u+1] = 0 ;
            for (i=1; i<=l[u+1]; i++)
              if ((s[u+1-i]!=-1) && (elp[u+1][i]!=0))
                d[u+1] ^= alpha_to_16[(s[u+1-i]+index_of_16[elp[u+1][i]])%nn] ;
            d[u+1] = index_of_16[d[u+1]] ;    /* put d[u+1] into index form */
          }
      } while ((u<nn-kk) && (l[u+1]<=tt)) ;

      u++ ;

      ////////////////////////////////////////////////////
      // error correction start!!!!!!!!!!!!!!!!!!!!!!!!!
      ////////////////////////////////////////////////////

      // CE 'or' SDC cases
      if (l[u]<=tt)         /* can correct error */
      {
/* put elp into index form */
         for (i=0; i<=l[u]; i++)   elp[u][i] = index_of_16[elp[u][i]] ;

/* find roots of the error location polynomial -> finding error location */
        for (i=1; i<=l[u]; i++)
          reg[i] = elp[u][i];
        
        count = 0 ;
        for (i=1; i<=nn; i++){  
          q = 1 ;
          for (j=1; j<=l[u]; j++)
            if (reg[j]!=-1){ 
              reg[j] = (reg[j]+j)%nn;
              q ^= alpha_to_16[reg[j]];
            }
          if (!q) {        /* store root and error location number indices */
            root[count] = i;
            loc[count] = nn-i; // -> error location
            count++ ;
          }
        }
        
        int CE_cases=1;
        //printf("error location check! (shortened RS code!)\n");
        for(int index=0; index<count; index++){
          if(loc[index]>=nn_short) // except for zero-padding part!!
            CE_cases=0;
        }
        

         // CE 'or' SDC cases
         if (count==l[u] && CE_cases==1){    /* no. roots = degree of elp hence <= tt errors*/
          //printf("CE 'or' SDC cases\n");
          //printf("count : %d\n",count);
/* form polynomial z(x) */
          for (i=1; i<=l[u]; i++){        /* Z[0] = 1 always - do not need */
            if ((s[i]!=-1) && (elp[u][i]!=-1))
              z[i] = alpha_to_16[s[i]] ^ alpha_to_16[elp[u][i]];
            else if ((s[i]!=-1) && (elp[u][i]==-1))
              z[i] = alpha_to_16[s[i]] ;
            else if ((s[i]==-1) && (elp[u][i]!=-1))
              z[i] = alpha_to_16[elp[u][i]] ;
            else
              z[i] = 0 ;
            for (j=1; j<i; j++){
              if ((s[j]!=-1) && (elp[u][i-j]!=-1))
                z[i] ^= alpha_to_16[(elp[u][i-j] + s[j])%nn];
            }
            z[i] = index_of_16[z[i]];         /* put into index form */
          }

  /* evaluate errors at locations given by error location numbers loc[i] */
          for (i=0; i<nn; i++){ 
            err[i] = 0;
            if (recd[i]!=-1)        /* convert recd[] to polynomial form */
              recd[i] = alpha_to_16[recd[i]];
            else  
              recd[i] = 0;
          }
          for (i=0; i<l[u]; i++){    /* compute numerator of error term first */
            err[loc[i]] = 1;       /* accounts for z[0] */
            for (j=1; j<=l[u]; j++){
              if (z[j]!=-1)
                err[loc[i]] ^= alpha_to_16[(z[j]+j*root[i])%nn];
            }
            if (err[loc[i]]!=0){
              err[loc[i]] = index_of_16[err[loc[i]]];
                q = 0;     /* form denominator of error term */
                for (j=0; j<l[u]; j++){
                  if (j!=i)
                    q += index_of_16[1^alpha_to_16[(loc[j]+root[i])%nn]];
                }
                q = q % nn;
                err[loc[i]] = alpha_to_16[(err[loc[i]]-q+nn)%nn];
                recd[loc[i]] ^= err[loc[i]];  /*recd[i] must be in polynomial form */
            }
          }

          return CE;
          /*
          printf("err (error values) : ");
          for(int index=0; index<nn; index++)
            printf("%d ",err[index]);
          printf("\n");
          */
        }
        // DUE cases
        else{    /* no. roots != degree of elp => >tt errors and cannot solve */
            for (i=0; i<nn; i++){ 
            err[i] = 0;
            if (recd[i]!=-1)        /* convert recd[] to polynomial form */
              recd[i] = alpha_to_16[recd[i]];
            else  
              recd[i] = 0;
          }
          return DUE;
        }
      }
    // DUE cases
     else{         /* elp has degree has degree >tt hence cannot solve */
         for (i=0; i<nn; i++){ 
            err[i] = 0;
            if (recd[i]!=-1)        /* convert recd[] to polynomial form */
              recd[i] = alpha_to_16[recd[i]];
            else  
              recd[i] = 0;
        }
       return DUE;
     }
  }
  // NE 'or' SDC cases
  else{       /* no non-zero syndromes => no errors: output received codeword */
    //printf("NE 'or' SDC cases!\n");
        for (i=0; i<nn; i++){ 
            err[i] = 0;
            if (recd[i]!=-1)        /* convert recd[] to polynomial form */
              recd[i] = alpha_to_16[recd[i]];
            else  
              recd[i] = 0;
        }
    return NE;
  }
}


int SDC_check(int BL, unsigned int Chip_array[][OECC_CW_LEN], int oecc_type, int recc_type)
{
    // Cacheline 에서 1이 남아있는지 검사
    // -> RECC 있으면 chip 0~9까지, RECC 없으면 chip 0~7까지

    int error_check=0;
    switch(recc_type){
        case Chipkill_correct_4_2:
            for(int Error_chip_pos=0; Error_chip_pos<CHIP_NUM; Error_chip_pos++){ // 0~9번째 chip까지
                for(int Error_pos=BL*4; Error_pos<(BL+2)*4; Error_pos++){
                    if(Chip_array[Error_chip_pos][Error_pos]==1){
                        error_check++;
                        return error_check;
                    }
                }
            }
            break;
        case Chipkill_correct_4_4:
            for(int Error_chip_pos=0; Error_chip_pos<CHIP_NUM; Error_chip_pos++){ // 0~9번째 chip까지
                for(int Error_pos=BL*4; Error_pos<(BL+4)*4; Error_pos++){
                    if(Chip_array[Error_chip_pos][Error_pos]==1){
                        error_check++;
                        return error_check;
                    }
                }
            }
            break;
        case Chipkill_correct_2_4:
            for(int Error_chip_pos=0; Error_chip_pos<CHIP_NUM; Error_chip_pos++){ // 0~9번째 chip까지
                for(int Error_pos=BL*4; Error_pos<(BL+4)*4; Error_pos++){ 
                    if(Chip_array[Error_chip_pos][Error_pos]==1){
                        error_check++;
                        return error_check;
                    }
                }
            }
            break;
        case Chipkill_correct_2_8:
            for(int Error_chip_pos=0; Error_chip_pos<CHIP_NUM; Error_chip_pos++){ // 0~9번째 chip까지
                for(int Error_pos=BL*4; Error_pos<(BL+8)*4; Error_pos++){ 
                    if(Chip_array[Error_chip_pos][Error_pos]==1){
                        error_check++;
                        return error_check;
                    }
                }
            }
            break;
        case Chipkill_correct_1_8:
            for(int Error_chip_pos=0; Error_chip_pos<CHIP_NUM; Error_chip_pos++){ // 0~9번째 chip까지
                for(int Error_pos=BL*4; Error_pos<(BL+8)*4; Error_pos++){
                    if(Chip_array[Error_chip_pos][Error_pos]==1){
                        error_check++;
                        return error_check;
                    }
                }
            }
            break;
        case Chipkill_correct_1_16:
            for(int Error_chip_pos=0; Error_chip_pos<CHIP_NUM; Error_chip_pos++){ // 0~9번째 chip까지
                for(int Error_pos=BL*4; Error_pos<(BL+16)*4; Error_pos++){
                    if(Chip_array[Error_chip_pos][Error_pos]==1){
                        error_check++;
                        return error_check;
                    }
                }
            }
            break;
        case NO_RLECC:
            for(int Error_chip_pos=0; Error_chip_pos<DATA_CHIP_NUM; Error_chip_pos++){ // 0~7번째 chip까지
                for(int Error_pos=0; Error_pos<OECC_DATA_LEN; Error_pos++){ // 0~127b 까지
                    if(Chip_array[Error_chip_pos][Error_pos]==1){
                        error_check++;
                        return error_check;
                    }
                }
            }
            break;
        default:
            break;
    }

    return error_check;
}


int main(int argc, char* argv[])
{
    // 1. generate the Galois field & compute the lookup table for this RS code
    generate_gf();

    // 2. H_Matrix 설정
    // SEC : OECC
    generator_oecc_H_matrix();

    // 3. 출력 파일 이름 설정 & oecc/fault/recc type 설정
    // main 함수의 argv parameter로 받는다.
    // run.py에서 변경 가능!!!

    // 파일명 예시 (OECC, RECC, Fault 순서)
    // ex : Bound_1_16_Chipkill_correct_4_2_SE_MBBE -> OECC에는 1 x 16 bounded SEC, RECC에는 4 x 2 Chipkill correct, Error는 2개 chip에서 각각 SE, MBBE가 발생한 경우
    // ex : Un_bound_Chipkill_correct_4_2_MBBE_SW -> OECC에는 un bounded SEC, RECC에는 4 x 2 Chipkill correct, Error는 2개 chip에서 각각 MBBE, SW가 발생한 경우

    string OECC="X", RECC="X", FAULT="X"; // => 파일 이름 생성을 위한 변수들. 그 이후로는 안쓰인다.
    int oecc_type, recc_type,fault_type; // => on-die ECC, Rank-level ECC, fault_type 분류를 위해 쓰이는 변수. 뒤에서도 계속 사용된다.
    oecc_recc_fault_type_assignment(OECC, FAULT, RECC, &oecc_type, &fault_type, &recc_type, atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
    
    string Result_file_name = OECC + "_" + RECC + "_" + FAULT + ".S";
    FILE *fp3=fopen(Result_file_name.c_str(),"w"); // c_str : string class에서 담고 있는 문자열을 c에서의 const char* 타입으로 변환하여 반환해주는 편리한 멤버함수

    // 4. 여기서부터 반복문 시작 (1억번)
    // DIMM 설정 (Channel에 있는 chip 구성을 기본으로 한다.)
    // DDR5 : x4 chip 기준으로 10개 chip이 있다. 각 chip은 on-die ECC codeword 136b이 있다.

    unsigned int Chip_array[CHIP_NUM][OECC_CW_LEN]; // 전체 chip 구성 (BL34 기준. [data : BL32, OECC-redundancy : BL2])
    int CE_cnt=0, DUE_cnt=0, SDC_cnt=0; // CE, DUE, SDC 횟수
    srand((unsigned int)time(NULL)); // 난수 시드값 계속 변화
    for(int runtime=0; runtime<RUN_NUM; runtime++){
        if(runtime%1000000==0){
            fprintf(fp3,"\n===============\n");
            fprintf(fp3,"Runtime : %d/%d\n",runtime,RUN_NUM);
            fprintf(fp3,"CE : %d\n",CE_cnt);
            fprintf(fp3,"DUE : %d\n",DUE_cnt);
            fprintf(fp3,"SDC : %d\n",SDC_cnt);
            fprintf(fp3,"\n===============\n");
	    fflush(fp3);
        }
        // 4-1. 10개 chip의 136b 전부를 0으로 초기화 (no-error)
        // 이렇게 하면 굳이 encoding을 안해도 된다. no-error라면 syndrome이 0으로 나오기 때문!
        for(int i=0; i<CHIP_NUM; i++)
            memset(Chip_array[i], 0, sizeof(unsigned int) * OECC_CW_LEN); 

        // 4-2. Error injection
        // [1] 2개의 chip을 선택 (Fault_Chip_position)
        // [2] error injection 각각 다른 위치의 chip 2개 뽑기 (0~9)
        // [3] chipkill은 각 bit마다 50% 확률로 bit/flip 발생
        vector<int> Fault_Chip_position;
        for (;;) {
            Fault_Chip_position.clear();
            Fault_Chip_position.push_back(rand()%CHIP_NUM); // Fault_Chip_position[0] // 0~9
            Fault_Chip_position.push_back(rand()%CHIP_NUM); // Fault_Chip_position[1] // 0~9
            if (Fault_Chip_position[0] != Fault_Chip_position[1]) break;
        }

        switch (fault_type){
            case SE_SE: // 1bit + 1 bit                
                error_injection_SE(Fault_Chip_position[0],Chip_array, oecc_type);
                error_injection_SE(Fault_Chip_position[1],Chip_array, oecc_type);
                break; 
            case SE_MBBE: // 1bit + (Multi bit Bounded Error injection) (More than 1 bit)
                error_injection_SE(Fault_Chip_position[0],Chip_array, oecc_type);
                error_injection_MBBE(Fault_Chip_position[1],Chip_array, oecc_type);
                break;
            case SE_SW: // 1bit + 1word
                error_injection_SE(Fault_Chip_position[0],Chip_array, oecc_type);
                error_injection_SW(Fault_Chip_position[1],Chip_array, oecc_type);
                break;
            case SE_SP: // 1bit + 1pin
                error_injection_SE(Fault_Chip_position[0],Chip_array, oecc_type);
                error_injection_SP(Fault_Chip_position[1],Chip_array, oecc_type);
                break;
            case SE_CHIPKILL: // 1bit + chipkill
                error_injection_SE(Fault_Chip_position[0],Chip_array, oecc_type);
                error_injection_CHIPKILL(Fault_Chip_position[1],Chip_array, oecc_type);
                break;
            case MBBE_MBBE: // (Multi bit Bounded Error injection) + (Multi bit Bounded Error injection) 
                error_injection_MBBE(Fault_Chip_position[0],Chip_array, oecc_type);
                error_injection_MBBE(Fault_Chip_position[1],Chip_array, oecc_type);
                break;
            case MBBE_SW: // (Multi bit Bounded Error injection) + 1word
                error_injection_MBBE(Fault_Chip_position[0],Chip_array, oecc_type);
                error_injection_SW(Fault_Chip_position[1],Chip_array, oecc_type);
                break;
            case MBBE_SP: // (Multi bit Bounded Error injection) + 1pin
                error_injection_MBBE(Fault_Chip_position[0],Chip_array, oecc_type);
                error_injection_SP(Fault_Chip_position[1],Chip_array, oecc_type);
                break;
            case MBBE_CHIPKILL: // (Multi bit Bounded Error injection) + chipkill
                error_injection_MBBE(Fault_Chip_position[0],Chip_array, oecc_type);
                error_injection_CHIPKILL(Fault_Chip_position[1],Chip_array, oecc_type);
                break;
            case SW_SW: // 1word + 1word
                error_injection_SW(Fault_Chip_position[0],Chip_array, oecc_type);
                error_injection_SW(Fault_Chip_position[1],Chip_array, oecc_type);
                break;
            case SW_SP: // 1word + 1pin
                error_injection_SW(Fault_Chip_position[0],Chip_array, oecc_type);
                error_injection_SP(Fault_Chip_position[1],Chip_array, oecc_type);
                break;
            case SW_CHIPKILL: // 1word + chipkill
                error_injection_SW(Fault_Chip_position[0],Chip_array, oecc_type);
                error_injection_CHIPKILL(Fault_Chip_position[1],Chip_array, oecc_type); 
                break;
            case SP_SP: // 1pin + 1pin
                error_injection_SP(Fault_Chip_position[0],Chip_array, oecc_type);
                error_injection_SP(Fault_Chip_position[1],Chip_array, oecc_type); 
                break;
            case SP_CHIPKILL: // 1pin + chipkill
                error_injection_SP(Fault_Chip_position[0],Chip_array, oecc_type);
                error_injection_CHIPKILL(Fault_Chip_position[1],Chip_array, oecc_type);
                break; 
            case CHIPKILL_CHIPKILL: // chipkill + chipkill
                error_injection_CHIPKILL(Fault_Chip_position[0],Chip_array, oecc_type);
                error_injection_CHIPKILL(Fault_Chip_position[1],Chip_array, oecc_type); 
                break;          
            default:
                break;
        }

        // 4-3. OECC
        // [1] Error를 넣은 chip에 대해서 SEC 실행
        // SEC : 136개의 1-bit error syndrome에 대응하면 correction 진행. 아니면 안함 (mis-correction을 최대한 막아보기 위함이다.)
        for(int Error_chip_pos=0; Error_chip_pos<CHIP_NUM; Error_chip_pos++) // 0~9
            error_correction_oecc(Error_chip_pos, Chip_array, oecc_type);

        /*
        4-4. RECC

        1. 각 BL2~8씩 묶어서 RECC 진행
            (1) Syndrome이 모두 0이면 => return NE
            (2) Syndrome이 0이 아니고, SSC의 syndrome과 일치하면 SSC 진행 [S1/S0이 a^0~a^9중 하나] => return CE
                -> 이때, error 발생 chip 위치는 CE 일때만 return 한다. (0~9). 이때만 의미가 있기 때문!
            (3) Syndrome이 0이 아니고, SSC의 syndrome과 일치하지 않으면 => return DUE
            (4) NE/CE return 했는데, 0이 아닌 값이 남아있으면 SDC (error 감지 못한 것이기 때문! 또는 mis-correction)

        2. Cacheline 단위로 (BL16) conservative 확인
            (1) CE : 8개 RECC 결과에서 SDC,DUE가 없고 전부 NE/CE 이고, error 발생 chip 위치(0~9)가 전부 같은 경우에만 CE로 처리 (NE는 고려 X. Error가 발생하지 않았기 때문!)
            (2) DUE : 8개 RECC 결과에서 SDC가 없고, 1개라도 DUE가 있는 경우 [나머지는 NE/CE]
                ** AMDCHIPKILL만 적용되는 조건! => 또는 전부 NE/CE이지만 error 발생 chip 위치가 다른 경우
            (3) SDC : 8개의 RECC 결과에서 SDC가 1개라도 있는 경우

        3. Cacheline 2개 (BL32) conservative 확인
            (1) CE : 2개 cacheline이 전부 CE인 경우
            (2) DUE : 2개 cacheline에서 SDC가 없고 1개라도 DUE가 있는 경우
            (3) SDC : 2개 cacheline에서 1개라도 SDC가 있는 경우
        
        */
        
        // 여기는 configuration 자동화 안 시켰음 (나중에 channel size, BL등을 바꿀거면 여기는 주의해서 보기!!)
        // set<int> error_chip_position; // RECC(BL2 단위)마다 error 발생한 chip 위치 저장 (CE인 경우에만 저장)
        int result_type_recc; // NE, CE, DUE, SDC 저장
        int final_result, final_result_1=CE,final_result_2=CE; // 각각 2개 cachline 고려한 최종 결과, 첫번째 cacheline, 두번째 cacheline 검사 결과
        int isConservative=0;
        switch(recc_type){
            case Chipkill_correct_4_2:{ // GF(2^8) 기준
                // 첫번째 cacheline
                for(int BL=0; BL<16; BL+=2){
                    int recd [255]={0,}, bb [2]={0,}; // codeword, redundancy
                    // recd에 codeworㅇ 옮기기
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            recd[index]^=(Chip_array[index][BL*4+symbol_index]<<symbol_index);
                    }

                    /* put recd[i] into index form */
                    for (int i=0; i<255; i++)
                        recd[i] = index_of_8[recd[i]];
                     // decoding
                    result_type_recc=decode_rs_8(recd, bb, 255, 253, 1, 10); //  recd, bb, nn, kk, tt, nn_short

                    // decoded codewrod를 Chip array에 update
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            Chip_array[index][BL*4+symbol_index]=recd[index]>>symbol_index;
                    }

                    // SDC 검사 (1이 남아있으면 SDC)
                    if(result_type_recc==CE || result_type_recc==NE){
                        int error_check=SDC_check(BL, Chip_array, oecc_type, recc_type);
                        if(error_check){
                            result_type_recc=SDC;
                        }
                    }
                    
                    // DUE 검사 (Restrained mode)
                    if(result_type_recc==DUE || final_result_1==DUE)
                        final_result_1=DUE;
                    else{ // 둘 중 우선순위가 큰 값 (SDC > CE > NE), 이전에 DUE가 나온 적이 없는 경우에만 들어갈 수 있다.
                        final_result_1 = (final_result_1>result_type_recc) ? final_result_1 : result_type_recc;
                    }
                }
                // 두번째 cacheline
                for(int BL=16; BL<32; BL+=2){ // BL : 16~31
                    int recd [255]={0,}, bb [2]={0,}; // codeword, redundancy
                    //encode_rs_8(recd, data, bb, 255, 253, 1); // parity 생성 (data가 all-zero인 경우), recd, data, bb, int nn, int kk, int tt
                    // recd에 codeword 옮기기
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            recd[index]^=(Chip_array[index][BL*4+symbol_index]<<symbol_index);
                    }
                    
                    /* put recd[i] into index form */
                    for (int i=0; i<255; i++)
                        recd[i] = index_of_8[recd[i]];
                     // decoding
                    result_type_recc=decode_rs_8(recd, bb, 255, 253, 1, 10); //  recd, bb, nn, kk, tt, nn_short

                    // decoded codewrod를 Chip array에 update
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            Chip_array[index][BL*4+symbol_index]=recd[index]>>symbol_index;
                    }

                    // SDC 검사 (1이 남아있으면 SDC)
                    if(result_type_recc==CE || result_type_recc==NE){
                        int error_check=SDC_check(BL, Chip_array, oecc_type, recc_type);
                        if(error_check){
                            result_type_recc=SDC;
                        }
                    }
                    
                    // DUE 검사 (Restrained mode)
                    if(result_type_recc==DUE || final_result_1==DUE)
                        final_result_1=DUE;
                    else{ // 둘 중 우선순위가 큰 값 (SDC > CE > NE), 이전에 DUE가 나온 적이 없는 경우에만 들어갈 수 있다.
                        final_result_1 = (final_result_1>result_type_recc) ? final_result_1 : result_type_recc;
                    }
                }

                if(final_result_2==NE || final_result_2==CE){
                    final_result_2 = (isConservative) ? DUE : CE;
                }


                // 2개 cacheline 비교해서 최종결과 update
                // SDC : 2개 cacheline 중에서 1개라도 SDC가 있으면 전체는 SDC
                // DUE : 2개 cacheline 중에서 SDC가 없고, 1개라도 DUE가 있으면 전체는 DUE
                // CE : 그 외 경우 (둘 다 CE)
                final_result = (final_result_1 > final_result_2) ? final_result_1 : final_result_2;
                }
                break;
            case Chipkill_correct_4_4:{ // GF(2^16) 기준
                // 첫번째 cacheline
                for(int BL=0; BL<16; BL+=4){
                    int recd [65535]={0,}, bb [2]={0,}; // codeword, redundancy
                    // recd에 codeword 옮기기
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            recd[index]^=(Chip_array[index][BL*4+symbol_index]<<symbol_index);
                    }

                    /* put recd[i] into index form */
                    for (int i=0; i<65535; i++)
                        recd[i] = index_of_16[recd[i]];
                    // decoding
                    result_type_recc=decode_rs_16(recd, bb, 65535, 65533, 1, 10); //  recd, bb, nn, kk, tt, nn_short

                    // decoded codewrod를 Chip array에 update
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            Chip_array[index][BL*4+symbol_index]=recd[index]>>symbol_index;
                    }

                    // SDC 검사 (1이 남아있으면 SDC)
                    if(result_type_recc==CE || result_type_recc==NE){
                        int error_check=SDC_check(BL, Chip_array, oecc_type, recc_type);
                        if(error_check){
                            result_type_recc=SDC;
                        }
                    }
                    
                    // DUE 검사 (Restrained mode)
                    if(result_type_recc==DUE || final_result_1==DUE)
                        final_result_1=DUE;
                    else{ // 둘 중 우선순위가 큰 값 (SDC > CE > NE), 이전에 DUE가 나온 적이 없는 경우에만 들어갈 수 있다.
                        final_result_1 = (final_result_1>result_type_recc) ? final_result_1 : result_type_recc;
                    }
                }
                // 두번째 cacheline
                for(int BL=16; BL<32; BL+=4){ // BL : 16~31
                    int recd [65535]={0,}, bb [2]={0,}; // codeword, redundancy
                    // recd에 codeword 옮기기
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            recd[index]^=(Chip_array[index][BL*4+symbol_index]<<symbol_index);
                    }
                    /* put recd[i] into index form */
                    for (int i=0; i<65535; i++)
                        recd[i] = index_of_16[recd[i]];
                    // decoding
                    result_type_recc=decode_rs_16(recd, bb, 65535, 65533, 1, 10);

                    // decoded codewrod를 Chip array에 update
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            Chip_array[index][BL*4+symbol_index]=recd[index]>>symbol_index;
                    }

                    // SDC 검사 (1이 남아있으면 SDC)
                    if(result_type_recc==CE || result_type_recc==NE){
                        int error_check=SDC_check(BL, Chip_array, oecc_type, recc_type);
                        if(error_check){
                            result_type_recc=SDC;
                        }
                    }
                    
                    // DUE 검사 (Restrained mode)
                    if(result_type_recc==DUE || final_result_1==DUE)
                        final_result_1=DUE;
                    else{ // 둘 중 우선순위가 큰 값 (SDC > CE > NE), 이전에 DUE가 나온 적이 없는 경우에만 들어갈 수 있다.
                        final_result_1 = (final_result_1>result_type_recc) ? final_result_1 : result_type_recc;
                    }
                }

                if(final_result_2==NE || final_result_2==CE){
                    final_result_2 = (isConservative) ? DUE : CE;
                }


                // 2개 cacheline 비교해서 최종결과 update
                // SDC : 2개 cacheline 중에서 1개라도 SDC가 있으면 전체는 SDC
                // DUE : 2개 cacheline 중에서 SDC가 없고, 1개라도 DUE가 있으면 전체는 DUE
                // CE : 그 외 경우 (둘 다 CE)
                final_result = (final_result_1 > final_result_2) ? final_result_1 : final_result_2;
                }
                break;
            case Chipkill_correct_2_4:{ // GF(2^8) 기준
                // 첫번째 cacheline
                for(int BL=0; BL<16; BL+=4){
                    int recd [255]={0,}, bb [4]={0,}; // codeword, redundancy
                    // recd에 codeword 옮기기
                    for(int index=0; index<10; index++){
                        recd[index]^=Chip_array[index][BL*4];
                        recd[index]^=Chip_array[index][BL*4+1]<<1;
                        recd[index]^=Chip_array[index][BL*4+4]<<2;
                        recd[index]^=Chip_array[index][BL*4+5]<<3;
                        recd[index]^=Chip_array[index][BL*4+8]<<4;
                        recd[index]^=Chip_array[index][BL*4+9]<<5;
                        recd[index]^=Chip_array[index][BL*4+12]<<6;
                        recd[index]^=Chip_array[index][BL*4+13]<<7;
                    }                    
                    for(int index=10; index<20; index++){
                        recd[index]^=Chip_array[index-10][BL*4+2];
                        recd[index]^=Chip_array[index-10][BL*4+3]<<1;
                        recd[index]^=Chip_array[index-10][BL*4+6]<<2;
                        recd[index]^=Chip_array[index-10][BL*4+7]<<3;
                        recd[index]^=Chip_array[index-10][BL*4+10]<<4;
                        recd[index]^=Chip_array[index-10][BL*4+11]<<5;
                        recd[index]^=Chip_array[index-10][BL*4+14]<<6;
                        recd[index]^=Chip_array[index-10][BL*4+15]<<7;
                    }

                    /* put recd[i] into index form */
                    for (int i=0; i<255; i++)
                        recd[i] = index_of_8[recd[i]];

                    // decoding
                    result_type_recc=decode_rs_8(recd, bb, 255, 251, 2, 20);

                    // decoded codewrod를 Chip array에 update
                    for(int index=0; index<10; index++){
                        Chip_array[index][BL*4]=recd[index];
                        Chip_array[index][BL*4+1]=recd[index]>>1;
                        Chip_array[index][BL*4+4]=recd[index]>>2;
                        Chip_array[index][BL*4+5]=recd[index]>>3;
                        Chip_array[index][BL*4+8]=recd[index]>>4;
                        Chip_array[index][BL*4+9]=recd[index]>>5;
                        Chip_array[index][BL*4+12]=recd[index]>>6;
                        Chip_array[index][BL*4+13]=recd[index]>>7;
                    }
                    for(int index=10; index<20; index++){
                        Chip_array[index-10][BL*4+2]=recd[index];
                        Chip_array[index-10][BL*4+3]=recd[index]>>1;
                        Chip_array[index-10][BL*4+6]=recd[index]>>2;
                        Chip_array[index-10][BL*4+7]=recd[index]>>3;
                        Chip_array[index-10][BL*4+10]=recd[index]>>4;
                        Chip_array[index-10][BL*4+11]=recd[index]>>5;
                        Chip_array[index-10][BL*4+14]=recd[index]>>6;
                        Chip_array[index-10][BL*4+15]=recd[index]>>7;
                    }

                    // SDC 검사 (1이 남아있으면 SDC)
                    if(result_type_recc==CE || result_type_recc==NE){
                        int error_check=SDC_check(BL, Chip_array, oecc_type, recc_type);
                        if(error_check){
                            result_type_recc=SDC;
                        }
                    }
                    
                    // DUE 검사 (Restrained mode)
                    if(result_type_recc==DUE || final_result_1==DUE)
                        final_result_1=DUE;
                    else{ // 둘 중 우선순위가 큰 값 (SDC > CE > NE), 이전에 DUE가 나온 적이 없는 경우에만 들어갈 수 있다.
                        final_result_1 = (final_result_1>result_type_recc) ? final_result_1 : result_type_recc;
                    }
                }
                // 두번째 cacheline
                for(int BL=16; BL<32; BL+=4){ // BL : 16~31
                    int recd [255]={0,}, bb [4]={0,}; // codeword, data, redundancy
                                        // recd에 codeword 옮기기
                    for(int index=0; index<10; index++){
                        recd[index]^=Chip_array[index][BL*4];
                        recd[index]^=Chip_array[index][BL*4+1]<<1;
                        recd[index]^=Chip_array[index][BL*4+4]<<2;
                        recd[index]^=Chip_array[index][BL*4+5]<<3;
                        recd[index]^=Chip_array[index][BL*4+8]<<4;
                        recd[index]^=Chip_array[index][BL*4+9]<<5;
                        recd[index]^=Chip_array[index][BL*4+12]<<6;
                        recd[index]^=Chip_array[index][BL*4+13]<<7;
                    }                    
                    for(int index=10; index<20; index++){
                        recd[index]^=Chip_array[index-10][BL*4+2];
                        recd[index]^=Chip_array[index-10][BL*4+3]<<1;
                        recd[index]^=Chip_array[index-10][BL*4+6]<<2;
                        recd[index]^=Chip_array[index-10][BL*4+7]<<3;
                        recd[index]^=Chip_array[index-10][BL*4+10]<<4;
                        recd[index]^=Chip_array[index-10][BL*4+11]<<5;
                        recd[index]^=Chip_array[index-10][BL*4+14]<<6;
                        recd[index]^=Chip_array[index-10][BL*4+15]<<7;
                    }

                    /* put recd[i] into index form */
                    for (int i=0; i<255; i++)
                        recd[i] = index_of_8[recd[i]];
                    // decoding
                    result_type_recc=decode_rs_8(recd, bb, 255, 251, 2, 20);

                    // decoded codewrod를 Chip array에 update
                    for(int index=0; index<10; index++){
                        Chip_array[index][BL*4]=recd[index];
                        Chip_array[index][BL*4+1]=recd[index]>>1;
                        Chip_array[index][BL*4+4]=recd[index]>>2;
                        Chip_array[index][BL*4+5]=recd[index]>>3;
                        Chip_array[index][BL*4+8]=recd[index]>>4;
                        Chip_array[index][BL*4+9]=recd[index]>>5;
                        Chip_array[index][BL*4+12]=recd[index]>>6;
                        Chip_array[index][BL*4+13]=recd[index]>>7;
                    }
                    for(int index=10; index<20; index++){
                        Chip_array[index-10][BL*4+2]=recd[index];
                        Chip_array[index-10][BL*4+3]=recd[index]>>1;
                        Chip_array[index-10][BL*4+6]=recd[index]>>2;
                        Chip_array[index-10][BL*4+7]=recd[index]>>3;
                        Chip_array[index-10][BL*4+10]=recd[index]>>4;
                        Chip_array[index-10][BL*4+11]=recd[index]>>5;
                        Chip_array[index-10][BL*4+14]=recd[index]>>6;
                        Chip_array[index-10][BL*4+15]=recd[index]>>7;
                    }

                    // SDC 검사 (1이 남아있으면 SDC)
                    if(result_type_recc==CE || result_type_recc==NE){
                        int error_check=SDC_check(BL, Chip_array, oecc_type, recc_type);
                        if(error_check){
                            result_type_recc=SDC;
                        }
                    }
                    
                    // DUE 검사 (Restrained mode)
                    if(result_type_recc==DUE || final_result_1==DUE)
                        final_result_1=DUE;
                    else{ // 둘 중 우선순위가 큰 값 (SDC > CE > NE), 이전에 DUE가 나온 적이 없는 경우에만 들어갈 수 있다.
                        final_result_1 = (final_result_1>result_type_recc) ? final_result_1 : result_type_recc;
                    }
                }

                if(final_result_2==NE || final_result_2==CE){
                    final_result_2 = (isConservative) ? DUE : CE;
                }


                // 2개 cacheline 비교해서 최종결과 update
                // SDC : 2개 cacheline 중에서 1개라도 SDC가 있으면 전체는 SDC
                // DUE : 2개 cacheline 중에서 SDC가 없고, 1개라도 DUE가 있으면 전체는 DUE
                // CE : 그 외 경우 (둘 다 CE)
                final_result = (final_result_1 > final_result_2) ? final_result_1 : final_result_2;
                }
                break;
            case Chipkill_correct_2_8:{ // GF(2^16) 기준
                // 첫번째 cacheline
                for(int BL=0; BL<16; BL+=8){
                    int recd [65535]={0,}, bb [4]={0,}; // codeword, data, redundancy
                    // recd에 codeword 옮기기
                    for(int index=0; index<10; index++){
                        recd[index]^=Chip_array[index][BL*4];
                        recd[index]^=Chip_array[index][BL*4+1]<<1;
                        recd[index]^=Chip_array[index][BL*4+4]<<2;
                        recd[index]^=Chip_array[index][BL*4+5]<<3;
                        recd[index]^=Chip_array[index][BL*4+8]<<4;
                        recd[index]^=Chip_array[index][BL*4+9]<<5;
                        recd[index]^=Chip_array[index][BL*4+12]<<6;
                        recd[index]^=Chip_array[index][BL*4+13]<<7;
                        recd[index]^=Chip_array[index][BL*4+16]<<8;
                        recd[index]^=Chip_array[index][BL*4+17]<<9;
                        recd[index]^=Chip_array[index][BL*4+20]<<10;
                        recd[index]^=Chip_array[index][BL*4+21]<<11;
                        recd[index]^=Chip_array[index][BL*4+24]<<12;
                        recd[index]^=Chip_array[index][BL*4+25]<<13;
                        recd[index]^=Chip_array[index][BL*4+28]<<14;
                        recd[index]^=Chip_array[index][BL*4+29]<<15;
                    }                    
                    for(int index=10; index<20; index++){
                        recd[index]^=Chip_array[index-10][BL*4+2];
                        recd[index]^=Chip_array[index-10][BL*4+3]<<1;
                        recd[index]^=Chip_array[index-10][BL*4+6]<<2;
                        recd[index]^=Chip_array[index-10][BL*4+7]<<3;
                        recd[index]^=Chip_array[index-10][BL*4+10]<<4;
                        recd[index]^=Chip_array[index-10][BL*4+11]<<5;
                        recd[index]^=Chip_array[index-10][BL*4+14]<<6;
                        recd[index]^=Chip_array[index-10][BL*4+15]<<7;
                        recd[index]^=Chip_array[index-10][BL*4+18]<<8;
                        recd[index]^=Chip_array[index-10][BL*4+19]<<9;
                        recd[index]^=Chip_array[index-10][BL*4+22]<<10;
                        recd[index]^=Chip_array[index-10][BL*4+23]<<11;
                        recd[index]^=Chip_array[index-10][BL*4+26]<<12;
                        recd[index]^=Chip_array[index-10][BL*4+27]<<13;
                        recd[index]^=Chip_array[index-10][BL*4+30]<<14;
                        recd[index]^=Chip_array[index-10][BL*4+31]<<15;
                    }

                    /* put recd[i] into index form */
                    for (int i=0; i<65535; i++)
                        recd[i] = index_of_16[recd[i]];
                    // decoding
                    result_type_recc=decode_rs_16(recd, bb, 65535, 65531, 2, 20);

                    // decoded codewrod를 Chip array에 update
                    for(int index=0; index<10; index++){
                        Chip_array[index][BL*4]=recd[index];
                        Chip_array[index][BL*4+1]=recd[index]>>1;
                        Chip_array[index][BL*4+4]=recd[index]>>2;
                        Chip_array[index][BL*4+5]=recd[index]>>3;
                        Chip_array[index][BL*4+8]=recd[index]>>4;
                        Chip_array[index][BL*4+9]=recd[index]>>5;
                        Chip_array[index][BL*4+12]=recd[index]>>6;
                        Chip_array[index][BL*4+13]=recd[index]>>7;
                        Chip_array[index][BL*4+16]=recd[index]>>8;
                        Chip_array[index][BL*4+17]=recd[index]>>9;
                        Chip_array[index][BL*4+20]=recd[index]>>10;
                        Chip_array[index][BL*4+21]=recd[index]>>11;
                        Chip_array[index][BL*4+24]=recd[index]>>12;
                        Chip_array[index][BL*4+25]=recd[index]>>13;
                        Chip_array[index][BL*4+28]=recd[index]>>14;
                        Chip_array[index][BL*4+29]=recd[index]>>15;
                    }
                    for(int index=10; index<20; index++){
                        Chip_array[index-10][BL*4+2]=recd[index];
                        Chip_array[index-10][BL*4+3]=recd[index]>>1;
                        Chip_array[index-10][BL*4+6]=recd[index]>>2;
                        Chip_array[index-10][BL*4+7]=recd[index]>>3;
                        Chip_array[index-10][BL*4+10]=recd[index]>>4;
                        Chip_array[index-10][BL*4+11]=recd[index]>>5;
                        Chip_array[index-10][BL*4+14]=recd[index]>>6;
                        Chip_array[index-10][BL*4+15]=recd[index]>>7;
                        Chip_array[index-10][BL*4+18]=recd[index]>>8;
                        Chip_array[index-10][BL*4+19]=recd[index]>>9;
                        Chip_array[index-10][BL*4+22]=recd[index]>>10;
                        Chip_array[index-10][BL*4+23]=recd[index]>>11;
                        Chip_array[index-10][BL*4+26]=recd[index]>>12;
                        Chip_array[index-10][BL*4+27]=recd[index]>>13;
                        Chip_array[index-10][BL*4+30]=recd[index]>>14;
                        Chip_array[index-10][BL*4+31]=recd[index]>>15;
                    }

                    // SDC 검사 (1이 남아있으면 SDC)
                    if(result_type_recc==CE || result_type_recc==NE){
                        int error_check=SDC_check(BL, Chip_array, oecc_type, recc_type);
                        if(error_check){
                            result_type_recc=SDC;
                        }
                    }
                    
                    // DUE 검사 (Restrained mode)
                    if(result_type_recc==DUE || final_result_1==DUE)
                        final_result_1=DUE;
                    else{ // 둘 중 우선순위가 큰 값 (SDC > CE > NE), 이전에 DUE가 나온 적이 없는 경우에만 들어갈 수 있다.
                        final_result_1 = (final_result_1>result_type_recc) ? final_result_1 : result_type_recc;
                    }
                }
                // 두번째 cacheline
                for(int BL=16; BL<32; BL+=8){ // BL : 16~31
                    int recd [65535]={0,}, bb [4]={0,}; // codeword, data, redundancy
                    // recd에 codeword 옮기기
                    for(int index=0; index<10; index++){
                        recd[index]^=Chip_array[index][BL*4];
                        recd[index]^=Chip_array[index][BL*4+1]<<1;
                        recd[index]^=Chip_array[index][BL*4+4]<<2;
                        recd[index]^=Chip_array[index][BL*4+5]<<3;
                        recd[index]^=Chip_array[index][BL*4+8]<<4;
                        recd[index]^=Chip_array[index][BL*4+9]<<5;
                        recd[index]^=Chip_array[index][BL*4+12]<<6;
                        recd[index]^=Chip_array[index][BL*4+13]<<7;
                        recd[index]^=Chip_array[index][BL*4+16]<<8;
                        recd[index]^=Chip_array[index][BL*4+17]<<9;
                        recd[index]^=Chip_array[index][BL*4+20]<<10;
                        recd[index]^=Chip_array[index][BL*4+21]<<11;
                        recd[index]^=Chip_array[index][BL*4+24]<<12;
                        recd[index]^=Chip_array[index][BL*4+25]<<13;
                        recd[index]^=Chip_array[index][BL*4+28]<<14;
                        recd[index]^=Chip_array[index][BL*4+29]<<15;
                    }                    
                    for(int index=10; index<20; index++){
                        recd[index]^=Chip_array[index-10][BL*4+2];
                        recd[index]^=Chip_array[index-10][BL*4+3]<<1;
                        recd[index]^=Chip_array[index-10][BL*4+6]<<2;
                        recd[index]^=Chip_array[index-10][BL*4+7]<<3;
                        recd[index]^=Chip_array[index-10][BL*4+10]<<4;
                        recd[index]^=Chip_array[index-10][BL*4+11]<<5;
                        recd[index]^=Chip_array[index-10][BL*4+14]<<6;
                        recd[index]^=Chip_array[index-10][BL*4+15]<<7;
                        recd[index]^=Chip_array[index-10][BL*4+18]<<8;
                        recd[index]^=Chip_array[index-10][BL*4+19]<<9;
                        recd[index]^=Chip_array[index-10][BL*4+22]<<10;
                        recd[index]^=Chip_array[index-10][BL*4+23]<<11;
                        recd[index]^=Chip_array[index-10][BL*4+26]<<12;
                        recd[index]^=Chip_array[index-10][BL*4+27]<<13;
                        recd[index]^=Chip_array[index-10][BL*4+30]<<14;
                        recd[index]^=Chip_array[index-10][BL*4+31]<<15;
                    }

                    /* put recd[i] into index form */
                    for (int i=0; i<65535; i++)
                        recd[i] = index_of_16[recd[i]];
                    // decoding
                    result_type_recc=decode_rs_16(recd, bb, 65535, 65531, 2, 20);

                    // decoded codewrod를 Chip array에 update
                    for(int index=0; index<10; index++){
                        Chip_array[index][BL*4]=recd[index];
                        Chip_array[index][BL*4+1]=recd[index]>>1;
                        Chip_array[index][BL*4+4]=recd[index]>>2;
                        Chip_array[index][BL*4+5]=recd[index]>>3;
                        Chip_array[index][BL*4+8]=recd[index]>>4;
                        Chip_array[index][BL*4+9]=recd[index]>>5;
                        Chip_array[index][BL*4+12]=recd[index]>>6;
                        Chip_array[index][BL*4+13]=recd[index]>>7;
                        Chip_array[index][BL*4+16]=recd[index]>>8;
                        Chip_array[index][BL*4+17]=recd[index]>>9;
                        Chip_array[index][BL*4+20]=recd[index]>>10;
                        Chip_array[index][BL*4+21]=recd[index]>>11;
                        Chip_array[index][BL*4+24]=recd[index]>>12;
                        Chip_array[index][BL*4+25]=recd[index]>>13;
                        Chip_array[index][BL*4+28]=recd[index]>>14;
                        Chip_array[index][BL*4+29]=recd[index]>>15;
                    }
                    for(int index=10; index<20; index++){
                        Chip_array[index-10][BL*4+2]=recd[index];
                        Chip_array[index-10][BL*4+3]=recd[index]>>1;
                        Chip_array[index-10][BL*4+6]=recd[index]>>2;
                        Chip_array[index-10][BL*4+7]=recd[index]>>3;
                        Chip_array[index-10][BL*4+10]=recd[index]>>4;
                        Chip_array[index-10][BL*4+11]=recd[index]>>5;
                        Chip_array[index-10][BL*4+14]=recd[index]>>6;
                        Chip_array[index-10][BL*4+15]=recd[index]>>7;
                        Chip_array[index-10][BL*4+18]=recd[index]>>8;
                        Chip_array[index-10][BL*4+19]=recd[index]>>9;
                        Chip_array[index-10][BL*4+22]=recd[index]>>10;
                        Chip_array[index-10][BL*4+23]=recd[index]>>11;
                        Chip_array[index-10][BL*4+26]=recd[index]>>12;
                        Chip_array[index-10][BL*4+27]=recd[index]>>13;
                        Chip_array[index-10][BL*4+30]=recd[index]>>14;
                        Chip_array[index-10][BL*4+31]=recd[index]>>15;
                    }

                    // SDC 검사 (1이 남아있으면 SDC)
                    if(result_type_recc==CE || result_type_recc==NE){
                        int error_check=SDC_check(BL, Chip_array, oecc_type, recc_type);
                        if(error_check){
                            result_type_recc=SDC;
                        }
                    }
                    
                    // DUE 검사 (Restrained mode)
                    if(result_type_recc==DUE || final_result_1==DUE)
                        final_result_1=DUE;
                    else{ // 둘 중 우선순위가 큰 값 (SDC > CE > NE), 이전에 DUE가 나온 적이 없는 경우에만 들어갈 수 있다.
                        final_result_1 = (final_result_1>result_type_recc) ? final_result_1 : result_type_recc;
                    }
                }

                if(final_result_2==NE || final_result_2==CE){
                    final_result_2 = (isConservative) ? DUE : CE;
                }


                // 2개 cacheline 비교해서 최종결과 update
                // SDC : 2개 cacheline 중에서 1개라도 SDC가 있으면 전체는 SDC
                // DUE : 2개 cacheline 중에서 SDC가 없고, 1개라도 DUE가 있으면 전체는 DUE
                // CE : 그 외 경우 (둘 다 CE)
                final_result = (final_result_1 > final_result_2) ? final_result_1 : final_result_2;
                }
                break;
            case Chipkill_correct_1_8:{ // GF(2^8) 기준
                // 첫번째 cacheline
                for(int BL=0; BL<16; BL+=8){
                    int recd [255]={0,}, bb [8]={0,}; // codeword, redundancy
                    // recd에 codeword 옮기기
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            recd[index]^=Chip_array[index][BL*4+symbol_index*4]<<symbol_index;
                    }
                    for(int index=10; index<20; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            recd[index]^=Chip_array[index-10][BL*4+symbol_index*4+1]<<symbol_index;
                    }
                    for(int index=20; index<30; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            recd[index]^=Chip_array[index-20][BL*4+symbol_index*4+2]<<symbol_index;
                    }
                    for(int index=30; index<40; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            recd[index]^=Chip_array[index-30][BL*4+symbol_index*4+3]<<symbol_index;
                    }

                    /* put recd[i] into index form */
                    for (int i=0; i<255; i++)
                        recd[i] = index_of_8[recd[i]];
                    // decoding
                    result_type_recc=decode_rs_8(recd, bb, 255, 247, 4, 40);

                    // decoded codewrod를 Chip array에 update
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            Chip_array[index][BL*4+symbol_index*4]=recd[index]>>symbol_index;
                    }
                    for(int index=10; index<20; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            Chip_array[index-10][BL*4+symbol_index*4+1]=recd[index]>>symbol_index;
                    }
                    for(int index=20; index<30; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            Chip_array[index-20][BL*4+symbol_index*4+2]=recd[index]>>symbol_index;
                    }
                    for(int index=30; index<40; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            Chip_array[index-30][BL*4+symbol_index*4+3]=recd[index]>>symbol_index;
                    }


                    // SDC 검사 (1이 남아있으면 SDC)
                    if(result_type_recc==CE || result_type_recc==NE){
                        int error_check=SDC_check(BL, Chip_array, oecc_type, recc_type);
                        if(error_check){
                            result_type_recc=SDC;
                        }
                    }
                    
                    // DUE 검사 (Restrained mode)
                    if(result_type_recc==DUE || final_result_1==DUE)
                        final_result_1=DUE;
                    else{ // 둘 중 우선순위가 큰 값 (SDC > CE > NE), 이전에 DUE가 나온 적이 없는 경우에만 들어갈 수 있다.
                        final_result_1 = (final_result_1>result_type_recc) ? final_result_1 : result_type_recc;
                    }
                }
                // 두번째 cacheline
                for(int BL=16; BL<32; BL+=8){ // BL : 16~31
                    int recd [255]={0,}, bb [8]={0,}; // codeword, redundancy
                    // recd에 codeword 옮기기
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            recd[index]^=Chip_array[index][BL*4+symbol_index*4]<<symbol_index;
                    }
                    for(int index=10; index<20; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            recd[index]^=Chip_array[index-10][BL*4+symbol_index*4+1]<<symbol_index;
                    }
                    for(int index=20; index<30; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            recd[index]^=Chip_array[index-20][BL*4+symbol_index*4+2]<<symbol_index;
                    }
                    for(int index=30; index<40; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            recd[index]^=Chip_array[index-30][BL*4+symbol_index*4+3]<<symbol_index;
                    }

                    /* put recd[i] into index form */
                    for (int i=0; i<255; i++)
                        recd[i] = index_of_8[recd[i]];
                    // decoding
                    result_type_recc=decode_rs_8(recd, bb, 255, 247, 4, 40);

                    // decoded codewrod를 Chip array에 update
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            Chip_array[index][BL*4+symbol_index*4]=recd[index]>>symbol_index;
                    }
                    for(int index=10; index<20; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            Chip_array[index-10][BL*4+symbol_index*4+1]=recd[index]>>symbol_index;
                    }
                    for(int index=20; index<30; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            Chip_array[index-20][BL*4+symbol_index*4+2]=recd[index]>>symbol_index;
                    }
                    for(int index=30; index<40; index++){
                        for(int symbol_index=0; symbol_index<8; symbol_index++)
                            Chip_array[index-30][BL*4+symbol_index*4+3]=recd[index]>>symbol_index;
                    }

                    // SDC 검사 (1이 남아있으면 SDC)
                    if(result_type_recc==CE || result_type_recc==NE){
                        int error_check=SDC_check(BL, Chip_array, oecc_type, recc_type);
                        if(error_check){
                            result_type_recc=SDC;
                        }
                    }
                    
                    // DUE 검사 (Restrained mode)
                    if(result_type_recc==DUE || final_result_1==DUE)
                        final_result_1=DUE;
                    else{ // 둘 중 우선순위가 큰 값 (SDC > CE > NE), 이전에 DUE가 나온 적이 없는 경우에만 들어갈 수 있다.
                        final_result_1 = (final_result_1>result_type_recc) ? final_result_1 : result_type_recc;
                    }
                }

                if(final_result_2==NE || final_result_2==CE){
                    final_result_2 = (isConservative) ? DUE : CE;
                }


                // 2개 cacheline 비교해서 최종결과 update
                // SDC : 2개 cacheline 중에서 1개라도 SDC가 있으면 전체는 SDC
                // DUE : 2개 cacheline 중에서 SDC가 없고, 1개라도 DUE가 있으면 전체는 DUE
                // CE : 그 외 경우 (둘 다 CE)
                final_result = (final_result_1 > final_result_2) ? final_result_1 : final_result_2;
                }
                break;
            case Chipkill_correct_1_16:{ // GF(2^16) 기준
                // 첫번째 cacheline
                for(int BL=0; BL<16; BL+=16){
                    int recd [65535]={0,}, bb [8]={0,}; // codeword, redundancy
                    // recd에 codeword 옮기기
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            recd[index]^=Chip_array[index][BL*4+symbol_index*4]<<symbol_index;
                    }
                    for(int index=10; index<20; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            recd[index]^=Chip_array[index-10][BL*4+symbol_index*4+1]<<symbol_index;
                    }
                    for(int index=20; index<30; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            recd[index]^=Chip_array[index-20][BL*4+symbol_index*4+2]<<symbol_index;
                    }
                    for(int index=30; index<40; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            recd[index]^=Chip_array[index-30][BL*4+symbol_index*4+3]<<symbol_index;
                    }

                    
                    /* put recd[i] into index form */
                    for (int i=0; i<65535; i++)
                        recd[i] = index_of_16[recd[i]];
                    // decoding
                    result_type_recc=decode_rs_16(recd, bb, 65535, 65527, 4, 40);

                    // decoded codewrod를 Chip array에 update
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            Chip_array[index][BL*4+symbol_index*4]=recd[index]>>symbol_index;
                    }
                    for(int index=10; index<20; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            Chip_array[index-10][BL*4+symbol_index*4+1]=recd[index]>>symbol_index;
                    }
                    for(int index=20; index<30; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            Chip_array[index-20][BL*4+symbol_index*4+2]=recd[index]>>symbol_index;
                    }
                    for(int index=30; index<40; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            Chip_array[index-30][BL*4+symbol_index*4+3]=recd[index]>>symbol_index;
                    }

                    // SDC 검사 (1이 남아있으면 SDC)
                    if(result_type_recc==CE || result_type_recc==NE){
                        int error_check=SDC_check(BL, Chip_array, oecc_type, recc_type);
                        if(error_check){
                            result_type_recc=SDC;
                        }
                    }
                    
                    // DUE 검사 (Restrained mode)
                    if(result_type_recc==DUE || final_result_1==DUE)
                        final_result_1=DUE;
                    else{ // 둘 중 우선순위가 큰 값 (SDC > CE > NE), 이전에 DUE가 나온 적이 없는 경우에만 들어갈 수 있다.
                        final_result_1 = (final_result_1>result_type_recc) ? final_result_1 : result_type_recc;
                    }
                }
                // 두번째 cacheline
                for(int BL=16; BL<32; BL+=16){ // BL : 16~31
                    int recd [65535]={0,}, bb [8]={0,}; // codeword, redundancy
                    // recd에 codeword 옮기기
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            recd[index]^=Chip_array[index][BL*4+symbol_index*4]<<symbol_index;
                    }
                    for(int index=10; index<20; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            recd[index]^=Chip_array[index-10][BL*4+symbol_index*4+1]<<symbol_index;
                    }
                    for(int index=20; index<30; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            recd[index]^=Chip_array[index-20][BL*4+symbol_index*4+2]<<symbol_index;
                    }
                    for(int index=30; index<40; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            recd[index]^=Chip_array[index-30][BL*4+symbol_index*4+3]<<symbol_index;
                    }

                    /* put recd[i] into index form */
                    for (int i=0; i<65535; i++)
                        recd[i] = index_of_16[recd[i]];
                    // decoding
                    result_type_recc=decode_rs_16(recd, bb, 65535, 65527, 4, 40);

                    // decoded codewrod를 Chip array에 update
                    for(int index=0; index<10; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            Chip_array[index][BL*4+symbol_index*4]=recd[index]>>symbol_index;
                    }
                    for(int index=10; index<20; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            Chip_array[index-10][BL*4+symbol_index*4+1]=recd[index]>>symbol_index;
                    }
                    for(int index=20; index<30; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            Chip_array[index-20][BL*4+symbol_index*4+2]=recd[index]>>symbol_index;
                    }
                    for(int index=30; index<40; index++){
                        for(int symbol_index=0; symbol_index<16; symbol_index++)
                            Chip_array[index-30][BL*4+symbol_index*4+3]=recd[index]>>symbol_index;
                    }

                    // SDC 검사 (1이 남아있으면 SDC)
                    if(result_type_recc==CE || result_type_recc==NE){
                        int error_check=SDC_check(BL, Chip_array, oecc_type, recc_type);
                        if(error_check){
                            result_type_recc=SDC;
                        }
                    }
                    
                    // DUE 검사 (Restrained mode)
                    if(result_type_recc==DUE || final_result_1==DUE)
                        final_result_1=DUE;
                    else{ // 둘 중 우선순위가 큰 값 (SDC > CE > NE), 이전에 DUE가 나온 적이 없는 경우에만 들어갈 수 있다.
                        final_result_1 = (final_result_1>result_type_recc) ? final_result_1 : result_type_recc;
                    }
                }

                if(final_result_2==NE || final_result_2==CE){
                    final_result_2 = (isConservative) ? DUE : CE;
                }


                // 2개 cacheline 비교해서 최종결과 update
                // SDC : 2개 cacheline 중에서 1개라도 SDC가 있으면 전체는 SDC
                // DUE : 2개 cacheline 중에서 SDC가 없고, 1개라도 DUE가 있으면 전체는 DUE
                // CE : 그 외 경우 (둘 다 CE)
                final_result = (final_result_1 > final_result_2) ? final_result_1 : final_result_2;
                }
                break;
            case NO_RLECC:{
                    int error_check = SDC_check(0, Chip_array, oecc_type, recc_type);
                    if(error_check) 
                        final_result=SDC;
                    else 
                        final_result=CE;
                }
                break;
            default:
                break;
        }

        // 4-5. CE/DUE/SDC 체크
        // 최종 update (2개 cacheline 전부 고려)
        // CE, DUE, SDC 개수 세기
        CE_cnt   += (final_result==CE)  ? 1 : 0;
        DUE_cnt  += (final_result==DUE) ? 1 : 0;
        SDC_cnt  += (final_result==SDC) ? 1 : 0;
            
    }
    // for문 끝!!

    // 최종 update
    fprintf(fp3,"\n===============\n");
    fprintf(fp3,"Runtime : %d\n",RUN_NUM);
    fprintf(fp3,"CE : %d\n",CE_cnt);
    fprintf(fp3,"DUE : %d\n",DUE_cnt);
    fprintf(fp3,"SDC : %d\n",SDC_cnt);
    fprintf(fp3,"\n===============\n");
    fflush(fp3);

    // 최종 update (소숫점 표현)
    fprintf(fp3,"\n===============\n");
    fprintf(fp3,"Runtime : %d\n",RUN_NUM);
    fprintf(fp3,"CE : %.11f\n",(double)CE_cnt/(double)RUN_NUM);
    fprintf(fp3,"DUE : %.11f\n",(double)DUE_cnt/(double)RUN_NUM);
    fprintf(fp3,"SDC : %.11f\n",(double)SDC_cnt/(double)RUN_NUM);
    fprintf(fp3,"\n===============\n");
    fflush(fp3);

    fclose(fp3);


    return 0;
}
