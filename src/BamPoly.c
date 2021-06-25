#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"

#include "BamCommonLibrary.h"




void	alignment_Poly(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t start, uint32_t end , uint32_t ref_length);

int	BamPoly(FILE *file_bam_i, FILE *file_bai_i,  FILE *file_bed_i, toolsFlags *ToolsFlags){

	uint8_t *stream_i;//[65536];
	uint8_t *stream_o;//[65536];
	uint8_t	*top;

	int	j;
	int	len_data;
	int	counter;
	int	ref_ID	= -1;


	z_stream infstream;
	bgfzTail BgfzTail;
	
	alignmentHeader	AlignmentHeader;

	uint8_t *address;


	uint64_t	offset_beg;
	uint64_t	offset_bgzf;
	uint64_t	offset_decomp;

	bamHeader	BamHeader;
	bedTable	BedTable;

	stream_i= calloc(65536,sizeof(uint8_t));
	stream_o= calloc(65536*2,sizeof(uint8_t));

	/* Initial*/
	stream_i[0] = 120;
	stream_i[1] = 156;
	/*---------------*/

	top = CatchBamHeader_v1 (file_bam_i, &BamHeader, stream_i, stream_o);
	counter	= top - stream_o;
	address = stream_o;

	/* comment 2020/03/30 YC
	if (ToolsFlags->flag_hide != 1){
		for (i = 0;i < BamHeader.n_ref;i++){
			printf("%s\t%d\n",BamHeader.chr_name[i], BamHeader.chr_length[i]);
		}
	}
	*/
	createRegionTable(file_bed_i, &BedTable, &BamHeader, ToolsFlags);
	
	printf("ID\tPATTERN\n");

	//BAI File
	for (j = 0; j < BamHeader.n_ref; j++){
		ref_ID = j;
		if (BedTable.table_end[ref_ID] != BedTable.table_start[ref_ID]){
			offset_beg = CatchOffset (file_bai_i, ref_ID, BedTable.start[BedTable.table_start[ref_ID]]);
			offset_bgzf = offset_beg >> 16;
			offset_decomp = offset_beg&65535;
			
			
			//Offset Bam File
			fseek(file_bam_i,offset_bgzf,SEEK_SET);
			len_data = BGFZBlock(file_bam_i);
			decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
			fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);

			address = stream_o + offset_decomp;
			counter = BgfzTail.I_size - offset_decomp;
			

			//Start Unzip and Process Bam File	
			while ( (len_data = BGFZBlock(file_bam_i)) > 0 ){
				memmove(stream_o, address, sizeof(uint8_t)*counter);

				decompressBlock(&infstream, stream_i, stream_o+counter, len_data, file_bam_i);
				fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);

				top	= stream_o + counter + BgfzTail.I_size;
				address = stream_o;

				if (top - address > (long int) sizeof(alignmentHeader)){
					memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
				}
				while ( top - address >= (AlignmentHeader.block_size + 4) ){
					if (AlignmentHeader.refID != ref_ID){
						break;
					}
					if ((uint32_t)AlignmentHeader.pos > BedTable.table_max[ref_ID]){
						//printf("Error\t%u\t%u\n", AlignmentHeader.pos, BedTable.start[BedTable.table_end[ref_ID]]);
						break;	
					}
					address += sizeof(alignmentHeader);
					
					//Core_START
					if ((AlignmentHeader.FLAG&1024)==0 || ToolsFlags->flag_dup==0 ){
						if ( ~ToolsFlags->flag_filter | SelectReads(address, &AlignmentHeader, 30, 30)){
							alignment_Poly(address, &AlignmentHeader, BedTable.start[0], BedTable.end[0], BamHeader.chr_length[ref_ID]);
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
				//printf("%u\t%u\t%u\n", BedTable.start[BedTable.table_start[ref_ID]], AlignmentHeader.pos, BedTable.table_max[ref_ID]);
				if (AlignmentHeader.refID != ref_ID){
					break;
				}
				if ((uint32_t) AlignmentHeader.pos > BedTable.table_max[ref_ID]){
					break;	
				}
				counter = top - address;
			}
			if (ToolsFlags->flag_hide == 1){
				printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1, BamHeader.n_ref, BamHeader.chr_name[ref_ID]);
			}
		}
	}
	
	free(stream_i);	
	free(stream_o);	
	return 0;
}

void	alignment_Poly(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t start, uint32_t end, uint32_t ref_length){

	int	i,j;

	uint32_t k;

	char	*read_name;
	uint32_t *cigar;
	uint8_t	*seq;

	int	op;
	int	op_len;
	char	residue;
	int	length;
	uint32_t	index;
	int	flag = 0;
//	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
	read_name	= (char *)stream;
	stream	+= AlignmentHeader->l_read_name;
//	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
	cigar	= (uint32_t *)stream;
	stream	+= ((AlignmentHeader->n_cigar_op)*sizeof(uint32_t));
//	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));
	seq	= (uint8_t *)stream;

	length = 0;
	// Maybe need to check pos >= 0 ?
	index = AlignmentHeader->pos;
	flag = 0;
	for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		if (CIGAR(op) == 'M'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			for (j = 0;j < op_len;j++){
				residue = Bin2SeqTop(seq[length >> 1],length&1);
				if (index >= start && index < end){
					if (flag == 0){
						printf("%s\t",read_name);	
						for (k=start;k <index;k++){
							printf(">");	
						}
					}
					printf("%c",residue);
					flag = 1;
				}
				length++;
				index++;
			}
		}else if (CIGAR(op) == 'I'){
			if (index >= start && index < end){
				if (flag == 0){
					printf("%s\t",read_name);	
					for (k=start;k <index;k++){
						printf(">");	
					}
				}
				printf("[");
				for (j = 0;j < op_len;j++){
					residue = Bin2SeqTop(seq[length >> 1],length&1);
					printf("%c",residue);
					flag = 1;
					length++;
				}
				printf("]");
			}else {
				length+=op_len;
			}
		}else if (CIGAR(op) == 'D'){
			for (j = 0;j < op_len;j++){
				if (index >= start && index < end){
					if (flag == 0){
						printf("%s\t",read_name);	
						for (k=start;k <index;k++){
							printf(">");	
						}
					}
					printf("-");
					flag = 1;
				}
				index++;
			}
		}else if (CIGAR(op) == 'S'){
			length += op_len;
		}else if (CIGAR(op) == 'N'){
			for (j = 0;j < op_len;j++){
				if (index >= start && index < end){
					if (flag == 0){
						printf("%s\t",read_name);	
						for (k=start;k <index;k++){
							printf(">");	
						}
					}
					printf("N");
					flag = 1;
				}
				index++;
			}
		//	index += op_len;
		}
		// CIGAR(op) == 'H' don't care about it.
		if (index >= end){
			break;
		}
	}
	if (flag == 1){
		for (k=index;k < end;k++){
			printf("<");	
		}
		printf("\n");	
	}
	
//	free(read_name);
//	free(cigar);
//	free(seq);
}
