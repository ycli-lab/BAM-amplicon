#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"


void	alignment_SinglePointQuality(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t *PosQuality, uint32_t position , uint32_t ref_length, uint8_t flag_mapq);

int	BamSinglePointQuality(FILE *file_bam_i, FILE *file_bai_i, toolsFlags *ToolsFlags){

	uint8_t *stream_i;//[65536];
	uint8_t *stream_o;//[65536];
	uint8_t	*top;

	int	i,j;
	int	len_data;
	int	counter;
	int	ref_ID	= -1;

	z_stream infstream;
	bgfzTail BgfzTail;
	
	alignmentHeader	AlignmentHeader;

	uint8_t *address;

	bamHeader	BamHeader;
	baiTable	BaiTable;
	
	uint64_t	offset_beg;
	uint64_t	offset_bgzf;
	uint64_t	offset_decomp;
	uint32_t	PosQuality[128*8];
	uint32_t	SumQuality[8];
	uint32_t	sumQuality = 0;
	
//	PosQuality = calloc(256,sizeof(uint32_t));
	memset(PosQuality,0,sizeof(uint32_t)*128*8);
	
	stream_i= calloc(65536,sizeof(uint8_t));
	stream_o= calloc(65536*2,sizeof(uint8_t));


	/* Initial*/
	stream_i[0] = 120;
	stream_i[1] = 156;
	/*---------------*/

	//Modify
		ToolsFlags->start--;	// (0-based)
	//

	top = CatchBamHeader_v1 (file_bam_i, &BamHeader, stream_i, stream_o);
	counter	= top - stream_o;
	address = stream_o;
	
	BuildBaiTable (file_bai_i, &BaiTable);


	for (i = 0;i < BamHeader.n_ref;i++){
		if (strcmp(ToolsFlags->chromosome, BamHeader.chr_name[i]) == 0){
			ref_ID = i;
//			printf("%s\t%d\n", BamHeader.chr_name[i], BamHeader.chr_length[i]);
			if( ToolsFlags->start + 1 > BamHeader.chr_length[i]){
				printf("[Error] Position is bigger than Size of %s\n",ToolsFlags->chromosome);
				return -1;	
			}
		}
	}

	/* comment 2020/03/30 YC	
	if (ToolsFlags->flag_hide == 1){
		for (i = 0;i < n_ref;i++){
			printf("%s\t%d\n",chr_name[i], chr_length[i]);
		}
	}
	*/


	//Offset Bam File
	offset_beg = ReturnOffset ( &BaiTable, ref_ID, ToolsFlags->start);
	offset_bgzf = offset_beg >> 16;
	offset_decomp = offset_beg&65535;

	if (offset_beg != 0){
		fseek(file_bam_i, offset_bgzf, SEEK_SET);	
		len_data = BGFZBlock(file_bam_i);
		decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
		
	
		fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
		top	= stream_o + BgfzTail.I_size;
		address = stream_o + offset_decomp;
		counter = BgfzTail.I_size - offset_decomp;
	}
	
	//Start Unzip and Process Bam File	
	while ( (len_data = BGFZBlock(file_bam_i)) > 0 ){
		memmove( stream_o, address, sizeof(uint8_t)*counter);
		decompressBlock(&infstream, stream_i, stream_o+counter, len_data, file_bam_i);
		fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);

		top	= stream_o + counter + BgfzTail.I_size;
		address = stream_o;

		if (top - address > (long int) sizeof(alignmentHeader)){
			memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
		}
		while ( top - address >= (AlignmentHeader.block_size + 4) ){

			if (AlignmentHeader.refID != ref_ID){
				if (ref_ID != -1){
					if (ToolsFlags->flag_hide == 1){
						printf("[Bam File Unzip %d / %d ] %s done\n", ref_ID+1, BamHeader.n_ref, BamHeader.chr_name[ref_ID]);
					}	
				}
				ref_ID = AlignmentHeader.refID;
			}

			if ((ref_ID == -1) || (strcmp(ToolsFlags->chromosome, BamHeader.chr_name[ref_ID]) != 0) || (AlignmentHeader.pos > ToolsFlags->start)){
				break;	
			}
	
			address += sizeof(alignmentHeader);
		
			//Core_START
			if ((AlignmentHeader.FLAG&1024)==0 || ToolsFlags->flag_dup==0 ){
				if (ToolsFlags->flag_mapq == 1){
					alignment_SinglePointQuality(address, &AlignmentHeader, PosQuality, ToolsFlags->start, BamHeader.chr_length[ref_ID], 1);
				}else {
					alignment_SinglePointQuality(address, &AlignmentHeader, PosQuality, ToolsFlags->start, BamHeader.chr_length[ref_ID], 0);
				}
			}
			//Core_END
			
			address += (AlignmentHeader.block_size + 4) - sizeof(alignmentHeader);

			if (top - address > (long int) sizeof(alignmentHeader)){
				memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
			}else {
				break;	
			}
		}
		counter = top - address;
		if ((ref_ID == -1) || (strcmp(ToolsFlags->chromosome, BamHeader.chr_name[ref_ID]) != 0) || (AlignmentHeader.pos > ToolsFlags->start )){
			break;	
		}
	}
	if (ref_ID != -1){
		if (ToolsFlags->flag_hide == 1){
			printf("[Bam File Unzip %d / %d ] %s done\n", ref_ID+1, BamHeader.n_ref, BamHeader.chr_name[ref_ID]);
		}
	}


	for (j = 0;j < 8;j++){
		SumQuality[j] = 0;
	}
	
	printf("Q.SCORE\tA_FOR\tC_FOR\tG_FOR\tT_FOR\tA_REV\tC_REV\tG_REV\tT_REV\t**(Q.SCORE: quality score; FOR: forward; REV: reverse)\n");
	for(i = 0;i < 64;i++){
		printf("%d\t",i);
		for (j = 0;j < 8;j++){
			printf("%u\t",PosQuality[i+j*128]);
			SumQuality[j] += i*PosQuality[i+j*128];
			sumQuality += i*PosQuality[i+j*128];
		}
		printf("\n");
	}
	printf("SUM\t");
	for (j = 0;j < 8;j++){
		printf("%u\t",SumQuality[j]);
	}
	printf("\n");




//	for(i = 128;i < 256;i++){
//		printf("%u\t",PosQuality[i]);
//		sumQuality += (i-128)*PosQuality[i];
//	}
//	printf("\n");
//	printf("%d\n", sumQuality);
	
	free(stream_i);	
	free(stream_o);	
	return 0;
}

void	alignment_SinglePointQuality(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t *PosQuality, uint32_t position, uint32_t ref_length, uint8_t flag_mapq){

	int	i;

	char	*read_name;
	uint32_t *cigar;
	uint8_t	*seq;
	char	*qual;

	uint32_t	op;
	uint32_t	op_len;
	char	residue;
	int	length;
	uint32_t	index;

	uint32_t	baseQuality;
	
	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));
	MemoryCopy(&stream,(void **)&qual	, AlignmentHeader->l_seq	, sizeof(char));

	length = 0;
	index = (uint32_t) AlignmentHeader->pos;
	
	for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		if (CIGAR(op) == 'M'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			if (index <= position && position < index+op_len){
				length += (position-index);
				residue = Bin2SeqTop(seq[length >> 1],length&1);
			
				//Care Mapping Qualuty
				if (flag_mapq == 1 && AlignmentHeader->MAPQ < qual[length]){
					baseQuality = AlignmentHeader->MAPQ;
				}else {
					baseQuality = qual[length];
				}
				//Care Strand
				if ((AlignmentHeader->FLAG&16)>>4){
					baseQuality += 512;
				}


				
				if ( residue == 'A'){
					PosQuality[baseQuality] += 1;
				}else if ( residue == 'C'){
					PosQuality[baseQuality+128] += 1;
				}else if ( residue == 'G'){
					PosQuality[baseQuality+256] += 1;
				}else if ( residue == 'T'){
					PosQuality[baseQuality+384] += 1;
				}else {

				}

				break;
			}
			index += op_len;
			length += op_len;
		}else if (CIGAR(op) == 'I'){
			length += op_len;
		}else if (CIGAR(op) == 'D'){
			index += op_len;
		}else if (CIGAR(op) == 'S'){
			length += op_len;
		}else if (CIGAR(op) == 'N'){
			index += op_len;
		}
		// CIGAR(op) == 'H' don't care about it.
		
	}
	
	free(read_name);
	free(cigar);
	free(seq);
	free(qual);
}
