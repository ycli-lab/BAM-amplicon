#ifndef BamCommon
#define BamCommon

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zlib.h"
#include <math.h>
#include "BamStruct.h"
#include <sys/time.h>




char	CIGAR(int A);
int	tag(void *stream);

void	BinToEncodedSeq_Top(uint8_t A,int MSB);
void	BinToEncodedSeq(uint8_t A);

void 	PrintAlignmenetHeader(alignmentHeader *AlignmentHeader);
void 	PrintSequence(uint8_t *seq,int32_t seq_length);
void 	PrintCIGAR(uint32_t *cigar,uint16_t n_cigar_op);
void 	PrintPosInformation(posInformation *PosInformation,int position);
int 	EndPosition(int32_t pos, uint16_t n_cigar_op, uint32_t *cigar);
void 	MemoryCopy(uint8_t **stream,void **ptr,int length,int size);
void 	PrintBlank(int size);
int 	compare_uint32_t (const void * a, const void * b);


void	strcat_triple(char *main, char *first, char *second, char *third, int length){
	memset(main,0,length);
	strcat(main,first);
	strcat(main,second);
	strcat(main,third);
}

void BinToEncodedSeq(uint8_t A){
	
	switch(A){
	case 0:	printf("=");	break;
	case 1:	printf("A");	break;
	case 2:	printf("C");	break;
	case 3:	printf("M");	break;
	case 4:	printf("G");	break;
	case 5:	printf("R");	break;
	case 6:	printf("S");	break;
	case 7:	printf("V");	break;
	case 8:	printf("T");	break;
	case 9:	printf("W");	break;
	case 10:printf("Y");	break;
	case 11:printf("H");	break;
	case 12:printf("K");	break;
	case 13:printf("D");	break;
	case 14:printf("B");	break;
	case 15:printf("N");	break;
	default:printf("Error: Binery to DNA\n");
		exit(1);
	}
}

void	BinToEncodedSeq_Top(uint8_t A,int MSB){
	uint8_t	low 	= A & 15;
	uint8_t	high	= (A & 240)>>4;
	if (MSB == 0){
		BinToEncodedSeq(high);
	} else {
		BinToEncodedSeq(low);
	}
}

void PrintSequence(uint8_t *seq,int32_t seq_length){
	int	i;
	for (i = 0;i < seq_length;i++){
		BinToEncodedSeq_Top(seq[i >> 1],i&1);
	}
}
void PrintCIGAR(uint32_t *cigar,uint16_t n_cigar_op){
	int	i;
	for (i = 0;i < n_cigar_op;i++){
		printf("%d%c",cigar[i] >> 4,CIGAR(cigar[i] & 7));
	}
}

int	tag(void *stream){
	char	tag[2];
	char	val_type;
	uint32_t value;

	memcpy(&tag,stream,2);
	memcpy(&val_type,stream+2,1);
	stream = stream + 3;
	
	switch(val_type){
	case 'A': printf("Unknown A\n");
		  exit(1);
		break;
	case 'c': memcpy(&value,stream,1);
		break;
	case 'C':  memcpy(&value,stream,1);
		break;
	case 's':  memcpy(&value,stream,2);
		break;
	case 'S':  memcpy(&value,stream,2);
		break;
	case 'i':  memcpy(&value,stream,4);
		break;
	case 'I':  memcpy(&value,stream,4);
		break;
	case 'f': printf("Unknown f\n");
		  exit(1);
		break;
	case 'Z': printf("Unknown Z\n");
		  exit(1);
		break;
	case 'H': printf("Unknown H\n");
		  exit(1);
		break;
	case 'B': printf("Unknown B\n");
		  exit(1);
		break;
	default: printf("Tag Error\n");
		exit(1);
	}
	
	printf("%c%c %c %d\n",tag[0],tag[1],val_type,value);
		
	return	(3+sizeof(value));
	
}

char	CIGAR(int A){
	switch(A){
	case 0:	return 'M';
	case 1:	return 'I';
	case 2:	return 'D';
	case 3:	return 'N';
	case 4:	return 'S';
	case 5:	return 'H';
	case 6:	return 'P';
	case 7:	return '=';
	case 8:	return 'X';
	default : printf("Error CIGAR\n");
		exit(1);
	}
}

void PrintAlignmenetHeader(alignmentHeader *AlignmentHeader){
//	printf("%d ",AlignmentHeader->block_size);
	printf("%2d ",AlignmentHeader->refID+1);
	printf("%10d ",AlignmentHeader->pos+1);
//	printf("%u ",AlignmentHeader->l_read_name);
//	printf("%u ",AlignmentHeader->MAPQ);
//	printf("%u ",AlignmentHeader->bin);
//	printf("%u ",AlignmentHeader->n_cigar_op);
	printf("%4u ",AlignmentHeader->FLAG);
	printf("%2d ",AlignmentHeader->l_seq);
//	printf("%d ",AlignmentHeader->next_refID+1);
//	printf("%d ",AlignmentHeader->nex_pos+1);
//	printf("%d ",AlignmentHeader->tlen);
}

void PrintPosInformation(posInformation *PosInformation,int position){
	printf("position:%10d ",position);
	printf("genotype:%c ",PosInformation->genotype);	
	printf("het:%c ",PosInformation->het);	
	printf("A:%5d ",PosInformation->A);	
	printf("T:%5d ",PosInformation->T);	
	printf("G:%5d ",PosInformation->G);	
	printf("C:%5d\n",PosInformation->C);	
}
void PrintBlank(int size){
	for (;size > 0;size--){
		printf(" ");
	}
}
int EndPosition(int32_t pos, uint16_t n_cigar_op, uint32_t *cigar){
	
	int	i;
	int	op;
	int	op_len;
//	int	flag = 1;
	int	end;

	end	= pos;
	for (i = 0;i < n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		if (CIGAR(op) == 'M'){
			end = end + op_len;	
		}else if (CIGAR(op) == 'I'){
		}else if (CIGAR(op) == 'D'){
			end = end + op_len;	
		}else if (CIGAR(op) == 'S'){
		}else if (CIGAR(op) == 'N'){
			end = end + op_len;	
		}
	}
	return (end-1);
}


void MemoryCopy(uint8_t **stream,void **ptr,int length,int size){
	*ptr	= calloc(length+1,size);
	memcpy(*ptr,*stream,length*size);
	*stream	+= length*size;
}
char Bin2Seq(int A){
	
	switch(A){
	case 0:	return '=';
	case 1:	return 'A';
	case 2: return 'C';
	case 3:	return 'M';
	case 4:	return 'G';
	case 5:	return 'R';
	case 6:	return 'S';
	case 7:	return 'V';
	case 8:	return 'T';
	case 9:	return 'W';
	case 10:return 'Y';
	case 11:return 'H';
	case 12:return 'K';
	case 13:return 'D';
	case 14:return 'B';
	case 15:return 'N';
	default:printf("Error: binery to DNA\n");
		exit(1);
	}
}

char Bin2SeqTop(uint8_t A,int MSB){
	int	low 	= A & 15;
	int	high	= (A & 240)>>4;
	if (MSB == 0){
		return Bin2Seq(high);
	} else{
		return Bin2Seq(low);
	}
}

void	decompressBlock(z_stream *infstream, uint8_t *stream_i, uint8_t *stream_o, int len_data, FILE *file_bam_i){

	int	err;
	infstream->zalloc 	= Z_NULL;
	infstream->zfree 	= Z_NULL;
	infstream->opaque 	= Z_NULL;
	infstream->avail_in 	= len_data+2;
	infstream->total_in 	= 0;
	infstream->next_in 	= stream_i;
	infstream->avail_out 	= 65536;
	infstream->total_out 	= 0;
	infstream->next_out 	= stream_o;

	fread(stream_i+2,sizeof(uint8_t),len_data,file_bam_i);
	

	err = inflateInit(infstream);
	if (err != Z_OK){	printf("inflateInit\n");	exit(1);	}
	err = inflate(infstream, Z_NO_FLUSH);
	if (err != Z_OK){	printf("inflate %d\n",err);	exit(1);	}
	err = inflateEnd(infstream);
	if (err != Z_OK){	printf("inflateEnd\n"); 	exit(1);	}

}



int BGFZBlock(FILE *file_ptr){

	bgfzHeader BGFZHeader;	

	if (fread(&BGFZHeader,sizeof(bgfzHeader),1,file_ptr) != 1){
		return -1;
	}
/*	printf("%d %d %d %d %d %d %d %d %d %d %d %d\n"
		,BGFZHeader.NUM_31
		,BGFZHeader.NUM_139
		,BGFZHeader.NUM_8
		,BGFZHeader.NUM_4
		,BGFZHeader.MTIME
		,BGFZHeader.XFL
		,BGFZHeader.OS
		,BGFZHeader.XLEN
		,BGFZHeader.NUM_66
		,BGFZHeader.NUM_67
		,BGFZHeader.NUM_2
		,BGFZHeader.BSIZE);
*/	
	return	((int) BGFZHeader.BSIZE - BGFZHeader.XLEN - 19);
}




int refInformation(void *stream, char *chr_name, int *chr_length){
	int32_t	l_name;
	char	*name;
	int32_t	l_ref;

	memcpy(&l_name,stream,sizeof(int32_t));	
	stream += 4;
	name = calloc(l_name+1,sizeof(char));
	memcpy(name,stream,sizeof(char)*l_name);
	stream += l_name;
	memcpy(&l_ref,stream,sizeof(int32_t));
	
	stream += 4;
	memcpy(chr_name,name,strlen(name));
	
	*chr_length = l_ref;
	free(name);
	return	(4+l_name+4);
}

void	alignment_Coverage(uint8_t *stream, alignmentHeader *AlignmentHeader, posCoverage *PosCoverage, uint32_t ref_length){

	int	i,j;

	char	*read_name;
	uint32_t *cigar;
	uint8_t	*seq;
//	char	*qual;

	int	op;
	int	op_len;
	char	residue;
	int	length;
	uint32_t	index;

	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));

	length = 0;
	index = AlignmentHeader->pos;
	for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		if (CIGAR(op) == 'M'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			for (j = 0;j < op_len;j++){
				residue = Bin2SeqTop(seq[length >> 1],length&1);

				if (residue == 'A'){		PosCoverage[index].A++;
				}else if (residue == 'T'){	PosCoverage[index].T++;
				}else if (residue == 'G'){	PosCoverage[index].G++;
				}else if (residue == 'C'){	PosCoverage[index].C++;
				}else if (residue == 'N'){	PosCoverage[index].N++;
				}else {
					printf("Error: Residue:%c\n",residue);
					exit(1);
				}
				length++;
				index++;
			}
		}else if (CIGAR(op) == 'I'){
			length += op_len;
		}else if (CIGAR(op) == 'D'){
			index += op_len;
		}else if (CIGAR(op) == 'S'){
			length += op_len;
		}else if (CIGAR(op) == 'N'){
			index += op_len;
		}
	}
	
	free(read_name);
	free(cigar);
	free(seq);
}

void	alignment_Quality(uint8_t *stream, alignmentHeader *AlignmentHeader, posQuality *PosQuality, uint32_t ref_length){

	int	i,j;

//	char	*read_name;
	uint32_t *cigar;
	uint8_t	*seq;
	char	*qual;
	uint8_t	base_qual;

	int	op;
	int	op_len;
//	char	residue;
	int	length;
	uint32_t	index;

//	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
	stream	+= AlignmentHeader->l_read_name;
	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));
	MemoryCopy(&stream,(void **)&qual	, AlignmentHeader->l_seq	, sizeof(char));
//	qual = calloc(AlignmentHeader->l_seq,sizeof(char));
///	memcpy(qual,stream,AlignmentHeader->l_seq);
//	for(i=0;i<AlignmentHeader->l_seq;i++){
//		qual[i]+=33;	
//	}

	length = 0;
	index = AlignmentHeader->pos;
	for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		if (CIGAR(op) == 'M'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			for (j = 0;j < op_len;j++){
//				residue = Bin2SeqTop(seq[length >> 1],length&1);
				base_qual = qual[length];
//				base_qual = qual[length] - 33;
//				if (index == 14368){
//					printf("%d\t%s\t%d\n",AlignmentHeader->MAPQ,read_name,AlignmentHeader->l_seq);
//					printf("%d,%d\n",base_qual,qual[length]);
//					getchar();
//				}
				PosQuality[index].cov_sum++;
				PosQuality[index].base_qual_sum+=(base_qual);
				PosQuality[index].base_qual_pow+=((base_qual)*(base_qual));
				PosQuality[index].map_qual_sum+=AlignmentHeader->MAPQ;
				PosQuality[index].map_qual_pow+=(AlignmentHeader->MAPQ*AlignmentHeader->MAPQ);
				length++;
				index++;
			}
		}else if (CIGAR(op) == 'I'){
			length += op_len;
		}else if (CIGAR(op) == 'D'){
			index += op_len;
		}else if (CIGAR(op) == 'S'){
			length += op_len;
		}else if (CIGAR(op) == 'N'){
			index += op_len;
		}
	}	
//	free(read_name);
	free(cigar);
	free(seq);
	free(qual);
}

void	alignment_Range(uint8_t *stream, alignmentHeader *AlignmentHeader, uint8_t *in_region, uint64_t *map_in, uint64_t *map_out, uint32_t ref_length){

	int	i,j;

//	char	*read_name;
	uint32_t *cigar;
//	uint8_t	*seq;
//	char	*qual;

	int	op;
	int	op_len;
	int	length;
	uint32_t	index;
	int	flag_map = 0;

	//MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
	stream	+= AlignmentHeader->l_read_name;
	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
	//MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));

	length = 0;
	index = AlignmentHeader->pos;
	for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		if (CIGAR(op) == 'M'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			for (j = 0;j < op_len;j++){
				if (in_region[index] == 1){
					*map_in = *map_in + 1;
					flag_map = 1;
					break;
				}
				length++;
				index++;
			}
		}else if (CIGAR(op) == 'I'){
			length += op_len;
		}else if (CIGAR(op) == 'D'){
			index += op_len;
		}else if (CIGAR(op) == 'S'){
			length += op_len;
		}else if (CIGAR(op) == 'N'){
			index += op_len;
		}

		if (flag_map == 1){
			break;	
		}
	}
	if (flag_map == 0){
		*map_out = *map_out + 1;	
	}
	
	//free(read_name);
	free(cigar);
	//free(seq);
}

void	alignment_List(uint8_t *stream, alignmentHeader *AlignmentHeader, FILE *file_o, char **chr_name){

	int	i;

	char	*read_name;
	uint32_t *cigar;
	uint8_t	*seq;
	char	*qual;

	int	op;
	int	op_len;
	char	residue;

	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));
	MemoryCopy(&stream,(void **)&qual	, AlignmentHeader->l_seq 	, sizeof(char));

	fprintf(file_o,"%s\t",read_name);
	fprintf(file_o,"%d\t",AlignmentHeader->FLAG);
	fprintf(file_o,"%s\t",chr_name[AlignmentHeader->refID]);
	fprintf(file_o,"%d\t",AlignmentHeader->pos+1);
	if (AlignmentHeader->next_refID >= 0){
		fprintf(file_o,"%s\t", chr_name[AlignmentHeader->next_refID]);
		fprintf(file_o,"%d\t", AlignmentHeader->next_pos+1);
	}else {	
		fprintf(file_o,"*\t");
		fprintf(file_o,"0\t");
	}
	fprintf(file_o,"%d\t",AlignmentHeader->tlen);

	for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		fprintf(file_o,"%d%c", op_len,CIGAR(op));
	}
	fprintf(file_o,"\t");
	for (i = 0;i < AlignmentHeader->l_seq;i++){
		residue = Bin2SeqTop(seq[i >> 1],i&1);
		fprintf(file_o,"%c", residue);
	}
	fprintf(file_o,"\t");
	fprintf(file_o,"\n");
	
	free(read_name);
	free(cigar);
	free(seq);
	free(qual);

}




//void PrintBaseDist (FILE *file_base, posCoverage *PosCoverage, uint8_t *in_region, char *chr_name, int chr_length, covBaseDistribution *BaseDist_All, int operation){
void PrintBaseDist (FILE *file_base, regionInformation *RegionInformation, bedTable *BedTable, char *chr_name, int chr_idx, covBaseDistribution *BaseDist_All, int operation){
	covBaseDistribution BaseDist;	
	uint32_t	i;
	int	j;
	int	j_start, j_end;
	uint32_t max_end	= 0;
	int	total;
	
	if (operation == 1){
		memset(&BaseDist, 0, sizeof(covBaseDistribution));
		
		for (i = BedTable->table_start[chr_idx]; i < BedTable->table_end[chr_idx]; i++){
			if (BedTable->start[i] > max_end){
				j_start	= 0;	// BedTable->start[i] - BedTable->start[i];
				j_end	= BedTable->end[i] - BedTable->start[i];
				max_end	= BedTable->end[i];
			}else if (BedTable->end[i] > max_end){
				j_start	= max_end - BedTable->start[i];
				j_end	= BedTable->end[i] - BedTable->start[i];
				max_end	= BedTable->end[i];
			}else {
				continue;
			}
			if (RegionInformation[i].count == 0){
				BaseDist.cov_0	+= (j_end-j_start);
				continue;
			}
			for(j = j_start; j < j_end; j++){
				total = RegionInformation[i].ALL[j];
				BaseDist.cov_all += total;
				if (total == 0){	BaseDist.cov_0++;
				}else if (total < 10){	BaseDist.cov_9++;
				}else if (total < 20){	BaseDist.cov_19++;
				}else if (total < 30){	BaseDist.cov_29++;
				}else if (total < 40){	BaseDist.cov_39++;
				}else if (total < 50){	BaseDist.cov_49++;
				}else if (total < 60){	BaseDist.cov_59++;
				}else if (total < 70){	BaseDist.cov_69++;
				}else if (total < 80){	BaseDist.cov_79++;
				}else if (total < 90){	BaseDist.cov_89++;
				}else if (total < 100){	BaseDist.cov_99++;
				}else {			BaseDist.cov_100++;
				}		
			}
		}
		BaseDist_All->cov_0	+=BaseDist.cov_0;
		BaseDist_All->cov_9	+=BaseDist.cov_9;
		BaseDist_All->cov_19	+=BaseDist.cov_19;
		BaseDist_All->cov_29	+=BaseDist.cov_29;
		BaseDist_All->cov_39	+=BaseDist.cov_39;
		BaseDist_All->cov_49	+=BaseDist.cov_49;
		BaseDist_All->cov_59	+=BaseDist.cov_59;
		BaseDist_All->cov_69	+=BaseDist.cov_69;
		BaseDist_All->cov_79	+=BaseDist.cov_79;
		BaseDist_All->cov_89	+=BaseDist.cov_89;
		BaseDist_All->cov_99	+=BaseDist.cov_99;
		BaseDist_All->cov_100	+=BaseDist.cov_100;
		BaseDist_All->cov_all	+=BaseDist.cov_all;
		fprintf(file_base,"%s\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n",
			chr_name,
			BaseDist.cov_0,
			BaseDist.cov_9,
			BaseDist.cov_19,
			BaseDist.cov_29,
			BaseDist.cov_39,
			BaseDist.cov_49,
			BaseDist.cov_59,
			BaseDist.cov_69,
			BaseDist.cov_79,
			BaseDist.cov_89,
			BaseDist.cov_99,
			BaseDist.cov_100);
	}else if (operation == 2){
		fprintf(file_base,"Total\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n",
			BaseDist_All->cov_0,
			BaseDist_All->cov_9,
			BaseDist_All->cov_19,
			BaseDist_All->cov_29,
			BaseDist_All->cov_39,
			BaseDist_All->cov_49,
			BaseDist_All->cov_59,
			BaseDist_All->cov_69,
			BaseDist_All->cov_79,
			BaseDist_All->cov_89,
			BaseDist_All->cov_99,
			BaseDist_All->cov_100,
			BaseDist_All->cov_all);
	}
}


//void PrintRegionDist (FILE *file_region_stat_o, FILE *file_region_ratio_o, bedTable *BedTable,  posCoverage *PosCoverage, bamHeader *BamHeader, int ref_ID, int threshold, covRegionDistribution *RegionDist_All, int operation){
void PrintRegionDist (FILE *file_region_stat_o, FILE *file_region_ratio_o, bedTable *BedTable, regionInformation *RegionInformation, bamHeader *BamHeader, int ref_ID, uint32_t threshold, covRegionDistribution *RegionDist_All, int operation){
	
	covRegionDistribution CovRegionDistribution;	
	uint32_t	i, j;
	uint32_t	length;

	uint32_t	*array_cov;
	uint32_t	index_cov;

//	on_target = calloc(BamHeader->chr_length[ref_ID], sizeof(uint8_t)); 
	
	if (operation == 1){
		memset(&CovRegionDistribution, 0,sizeof(covRegionDistribution));


		//Calculate the number of target region

		//Calculate the ave coverage
		for (i = BedTable->table_start[ref_ID]; i < BedTable->table_end[ref_ID]; i++){
			
			if ( RegionInformation[i].count == 0){
				CovRegionDistribution.per_0++;
				RegionInformation[i].valid_percnt	= 0.0;
				RegionInformation[i].ave_cov	= 0;
				continue;
			}
			
			length = RegionInformation[i].length;
			RegionInformation[i].all_base	= 0;
			RegionInformation[i].valid_length	= 0;	
			for (j = 0; j < length; j++){
//				on_target[j+BedTable->start[i]]++;
				if ( RegionInformation[i].ALL[j] > threshold){
					RegionInformation[i].all_base += RegionInformation[i].ALL[j];
					RegionInformation[i].valid_length++;
				}
			}
			RegionInformation[i].valid_percnt	= ((double)RegionInformation[i].valid_length)/length;

			CovRegionDistribution.per_all += RegionInformation[i].valid_percnt;
			if (RegionInformation[i].valid_percnt == 0){		CovRegionDistribution.per_0++;	
			}else if (RegionInformation[i].valid_percnt < 0.5){	CovRegionDistribution.per_50++;	
			}else if (RegionInformation[i].valid_percnt < 0.9){	CovRegionDistribution.per_90++;	
			}else if (RegionInformation[i].valid_percnt != length){	CovRegionDistribution.per_99++;	
			}else {				CovRegionDistribution.per_100++;	
			}
			
			if (RegionInformation[i].valid_length == 0){
				RegionInformation[i].ave_cov = 0.0;
			}else {
				RegionInformation[i].ave_cov = ((float)RegionInformation[i].all_base)/RegionInformation[i].valid_length;	
			}
		}
		

		for (i = BedTable->table_start[ref_ID]; i < BedTable->table_end[ref_ID]; i++){
			if ( RegionInformation[i].count == 0){
				RegionInformation[i].cov_Q3	= 0;
				RegionInformation[i].cov_Q2	= 0;
				continue;
			}
			length	= RegionInformation[i].length;
			array_cov = calloc( length, sizeof(uint32_t));
			index_cov = 0;
			for (j = 0 ;j < RegionInformation[i].length; j++){
				if (RegionInformation[i].overlap[j] == 0){
//				if (on_target[j+BedTable->start[i]] == 1){
					array_cov[index_cov]	= RegionInformation[i].ALL[j];
					index_cov++;
				}else {
					length--;
				}
			}
			qsort ( array_cov, length, sizeof(uint32_t), compare_uint32_t);
			RegionInformation[i].cov_Q3	= array_cov[(int)(0.5 + (float)length * 0.75)];
			RegionInformation[i].cov_Q2	= array_cov[(int)(0.5 + (float)length * 0.5)];
			free(array_cov);
		}



		for (i = BedTable->table_start[ref_ID]; i < BedTable->table_end[ref_ID]; i++){
			fprintf(file_region_ratio_o,"%s\t%u\t%u\t%u\t%u\t%u\t%u\t%f\t%f\t%u\t%u",
				BamHeader->chr_name[ref_ID],
				BedTable->start[i],
				BedTable->end[i],
				RegionInformation[i].length,
				RegionInformation[i].count,
				RegionInformation[i].length-RegionInformation[i].valid_length,  
				RegionInformation[i].all_base,
				RegionInformation[i].ave_cov, 
				RegionInformation[i].valid_percnt,
				RegionInformation[i].cov_Q2,
				RegionInformation[i].cov_Q3);
			fprintf(file_region_ratio_o,"\n");
		}
		
		RegionDist_All->per_0 += CovRegionDistribution.per_0;
		RegionDist_All->per_50 += CovRegionDistribution.per_50;
		RegionDist_All->per_90 += CovRegionDistribution.per_90;
		RegionDist_All->per_99 += CovRegionDistribution.per_99;
		RegionDist_All->per_100 += CovRegionDistribution.per_100;
		RegionDist_All->per_all += CovRegionDistribution.per_all;
	//	printf("RegionDist_All->per_all:\t%lf\n",RegionDist_All->per_all);
		fprintf(file_region_stat_o, "%s\t%d\t%d\t%d\t%d\t%d\n",
			BamHeader->chr_name[ref_ID],
			CovRegionDistribution.per_0,
			CovRegionDistribution.per_50,
			CovRegionDistribution.per_90,
			CovRegionDistribution.per_99,
			CovRegionDistribution.per_100);



	}else if (operation == 2){
		fprintf(file_region_stat_o, "Total\t%d\t%d\t%d\t%d\t%d\t%lf\n",
			RegionDist_All->per_0,
			RegionDist_All->per_50,
			RegionDist_All->per_90,
			RegionDist_All->per_99,
			RegionDist_All->per_100,
			RegionDist_All->per_all);
	}
}


void	alignment_DepthTxt(uint8_t *stream, alignmentHeader *AlignmentHeader, posCoverage *PosCoverage, uint32_t ref_length){

	int	i,j;

	uint32_t *cigar;

	int	op;
	int	op_len;
	int	length;
	uint32_t	index;


	//MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
	stream	+= AlignmentHeader->l_read_name;
	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
	//MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));

	length = 0;
	index = AlignmentHeader->pos;

	for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		if (CIGAR(op) == 'M'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			for (j = 0;j < op_len;j++){
				PosCoverage[index].ALL++;;
				index++;
			}
			length += op_len;
		}else if (CIGAR(op) == 'I'){
			length += op_len;
		}else if (CIGAR(op) == 'D'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			for (j = 0;j < op_len;j++){
//				PosCoverage[index].DEL++;;
				PosCoverage[index].ALL++;;
				index++;
			}
		}else if (CIGAR(op) == 'S'){
			length += op_len;
		}else if (CIGAR(op) == 'N'){
			index += op_len;
		}
		// CIGAR(op) == 'H' don't care about it.
		
	}
	
	//free(read_name);
	free(cigar);
	//free(seq);

}

void	alignment_CoverageTxt(uint8_t *stream, alignmentHeader *AlignmentHeader, posCoverage *PosCoverage, uint32_t *Coverage, uint32_t ref_length){

	int	i,j;

	char	*read_name;
	uint32_t *cigar;
	uint8_t	*seq;

	int	op;
	int	op_len;
	char	residue;
	int	length;
	uint32_t	index;


	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));

	length = 0;
	index = AlignmentHeader->pos;

	for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		if (CIGAR(op) == 'M'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			for (j = 0;j < op_len;j++){
				residue = Bin2SeqTop(seq[length >> 1],length&1);
				if (residue == 'A'){		PosCoverage[index].A++;
				}else if (residue == 'T'){	PosCoverage[index].T++;
				}else if (residue == 'G'){	PosCoverage[index].G++;
				}else if (residue == 'C'){	PosCoverage[index].C++;
				}else if (residue == 'N'){	PosCoverage[index].N++;
				}else {
					printf("Error: Residue:%c\n",residue);
					exit(1);
				}
				Coverage[index]++;
				length++;
				index++;
			}
		}else if (CIGAR(op) == 'I'){
			length += op_len;
		}else if (CIGAR(op) == 'D'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			for (j = 0;j < op_len;j++){
				Coverage[index]++;
				index++;
			}
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

}

uint32_t	alignment_EndPosition(uint8_t *stream, alignmentHeader *AlignmentHeader, uint32_t ref_length){

	int	i;

	uint32_t *cigar;

	int	op;
	int32_t	op_len;
	int	length;
	uint32_t	index;


//	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
//	printf("%.*s\n", AlignmentHeader->l_read_name, stream);
	stream	+= AlignmentHeader->l_read_name;
	
//	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
	cigar	= (uint32_t *)stream;
	stream	+= AlignmentHeader->n_cigar_op*sizeof(uint32_t);
//	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));

	length = 0;
	index = AlignmentHeader->pos;

	for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		if (CIGAR(op) == 'M'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			length	+= op_len;
			index	+= op_len;
		}else if (CIGAR(op) == 'I'){
			length += op_len;
		}else if (CIGAR(op) == 'D'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			index += op_len;
		}else if (CIGAR(op) == 'S'){
			length += op_len;
		}else if (CIGAR(op) == 'N'){
			index += op_len;
		}
		// CIGAR(op) == 'H' don't care about it.
		
	}
	//free(read_name);
//	free(cigar);
//	free(seq);
	return index;	

}
int	alignment_DepthDist_TR(uint8_t *stream, alignmentHeader *AlignmentHeader, posCoverage *PosCoverage, uint32_t start, uint32_t end, uint32_t pos_offset){
//void	alignment_DepthDist_TR(uint8_t *stream, alignmentHeader *AlignmentHeader, posCoverage *PosCoverage, uint32_t start, uint32_t end){

	int	i,j;

//	char	*read_name;
	uint32_t *cigar;
	uint8_t	*seq;
//	char	*qual;

	int	op;
	int32_t	op_len;
	int32_t	_op_len;
	char	residue;
	int	length;
	uint32_t	index;


//	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
	stream	+= AlignmentHeader->l_read_name;
	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));

	length = 0;
	index = (uint32_t) AlignmentHeader->pos;

	for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		_op_len	= op_len;
		if (CIGAR(op) == 'M'){
			if (index < start){
//				if (op_len > (start - index)){
				if (index + op_len > start){
					length	+= (start - index);
					op_len	-= (start - index);
					_op_len	-= (start - index);
					index	= start;
				}else {
					length	+= op_len;
					index	+= op_len;
					op_len	= 0;
					_op_len	= 0;
				}
			}
			if (index > end){
				op_len	= 0;
			}else if (index + op_len > end){
				op_len = end - index;
			}
			for (j = 0;j < op_len;j++){
				residue = Bin2SeqTop(seq[length >> 1],length&1);
				if (residue == 'A'){		PosCoverage[index-pos_offset].A++;
				}else if (residue == 'T'){	PosCoverage[index-pos_offset].T++;
				}else if (residue == 'G'){	PosCoverage[index-pos_offset].G++;
				}else if (residue == 'C'){	PosCoverage[index-pos_offset].C++;
				}else if (residue == 'N'){	PosCoverage[index-pos_offset].N++;
				}else {
					printf("Error: Residue:%c\n",residue);
					exit(1);
				}
				length++;
				index++;
			}
			index	+= (_op_len - op_len);
		}else if (CIGAR(op) == 'I'){
			length += op_len;
		}else if (CIGAR(op) == 'D'){
			if (index < start){
				if (index + op_len > start){
//				if (op_len > (start - index)){
					op_len	-= (start - index);
					_op_len	-= (start - index);
					index	= start;
				}else {
					index	+= op_len;
					op_len	= 0;
					_op_len	= 0;
				}
			}
			if (index > end){
				op_len	= 0;
			}else if (index + op_len > end){
				op_len = end - index;
			}
			for (j = 0;j < op_len;j++){
				PosCoverage[index-pos_offset].DEL++;
				index++;
			}
			index	+= (_op_len - op_len);
		}else if (CIGAR(op) == 'S'){
			length += op_len;
		}else if (CIGAR(op) == 'N'){
			index += op_len;
		}
		// CIGAR(op) == 'H' don't care about it.
		
	}
	
	//free(read_name);
	free(cigar);
	free(seq);
	return index;
}

void	alignment_DepthDist(uint8_t *stream, alignmentHeader *AlignmentHeader, posCoverage *PosCoverage, uint32_t ref_length){

	int	i,j;

//	char	*read_name;
	uint32_t *cigar;
	uint8_t	*seq;
//	char	*qual;

	int	op;
	int	op_len;
	char	residue;
	int	length;
	uint32_t	index;


//	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
	stream	+= AlignmentHeader->l_read_name;
//	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
	cigar	= (uint32_t *)stream;
	stream	+= AlignmentHeader->n_cigar_op*sizeof(uint32_t);
//	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));
	seq	= (uint8_t *)stream;

	length = 0;
	index = (uint32_t) AlignmentHeader->pos;

	for (i = 0;i < AlignmentHeader->n_cigar_op;i++){
		op	= cigar[i] & 7;
		op_len	= cigar[i] >> 4;
		if (CIGAR(op) == 'M'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			for (j = 0;j < op_len;j++){
				residue = Bin2SeqTop(seq[length >> 1],length&1);
				if (residue == 'A'){		PosCoverage[index].A++;
				}else if (residue == 'T'){	PosCoverage[index].T++;
				}else if (residue == 'G'){	PosCoverage[index].G++;
				}else if (residue == 'C'){	PosCoverage[index].C++;
				}else if (residue == 'N'){	PosCoverage[index].N++;
				}else {
					printf("Error: Residue:%c\n",residue);
					exit(1);
				}
				length++;
				index++;
			}
		}else if (CIGAR(op) == 'I'){
			length += op_len;
		}else if (CIGAR(op) == 'D'){
			if (index + op_len > ref_length){
				op_len = ref_length - index;
			}
			for (j = 0;j < op_len;j++){
				PosCoverage[index].DEL++;
				index++;
			}
		}else if (CIGAR(op) == 'S'){
			length += op_len;
		}else if (CIGAR(op) == 'N'){
			index += op_len;
		}
		// CIGAR(op) == 'H' don't care about it.
		
	}
	
	//free(read_name);
	//free(cigar);
	//free(seq);
}

int compare_uint32_t (const void * a, const void * b)
{
	  return ( *(uint32_t*)a - *(uint32_t*)b );
}

void BuildBaiTable (FILE *file_bai_i, baiTable *BaiTable ){

	uint32_t	i,j,m;
	char	BAMI[4];
	uint32_t	n_bin;
	uint32_t	bin;
	uint32_t	n_chunk;
	uint64_t	chunk_beg;
	uint64_t	chunk_end;
	
	//BAI File
	
	rewind(file_bai_i);
	fread(BAMI,4, sizeof(char), file_bai_i);

	fread(&BaiTable->n_ref_bai, 1, sizeof(int32_t), file_bai_i);

	BaiTable->ioffset	= calloc(sizeof(uint64_t *), BaiTable->n_ref_bai);
	BaiTable->n_intv = calloc(sizeof(int32_t), BaiTable->n_ref_bai);

	for(i = 0;i < BaiTable->n_ref_bai;i++){

		fread(&n_bin,1,sizeof(int32_t),file_bai_i);
		for (j = 0;j < n_bin;j++){
			fread(&bin,1,sizeof(uint32_t),file_bai_i);
			fread(&n_chunk,1,sizeof(int32_t),file_bai_i);
//			fseek(file_bai_i, n_chunk<<4, SEEK_CUR);
//			printf("%d\n", n_chunk);
			for(m = 0;m < n_chunk;m++){
				fread(&chunk_beg,1,sizeof(uint64_t),file_bai_i);
				fread(&chunk_end,1,sizeof(uint64_t),file_bai_i);
			}
		}

		fread(&BaiTable->n_intv[i], 1, sizeof(int32_t), file_bai_i);
		BaiTable->ioffset[i] = calloc(sizeof(uint64_t), BaiTable->n_intv[i]);
//		printf("%d\n", BaiTable->n_intv[i]);
		fread(BaiTable->ioffset[i], BaiTable->n_intv[i], sizeof(uint64_t), file_bai_i);
	}
}

uint64_t ReturnOffset (baiTable *BaiTable, int chr_ID, uint32_t start ){

	int	i;
	
	uint64_t	offset_beg;
//	uint64_t	offset_bgzf;
//	uint64_t	offset_decomp;
	uint32_t	start_index;

	//BAI File

	i = chr_ID;
	start_index	= start >> 14;
	if (start_index >= BaiTable->n_intv[i]){
//		printf("[Warning][ReturnOffset] start_index (%u) over BaiTable->n_intv[%u] (%u) \n", start_index, i, BaiTable->n_intv[i]);
		return 0;
	}

	offset_beg = BaiTable->ioffset[i][start_index];

	while ( offset_beg == 0 && (chr_ID != 0 || start_index > 0)){
		if (start_index == 0){
			i--;
			if (i < 0 || BaiTable->n_intv[i] == 0){
				return 0;
			}
			start_index = BaiTable->n_intv[i] - 1;
		}else {
			start_index --;
		}
		offset_beg = BaiTable->ioffset[i][start_index];
	}
	return offset_beg;
}

uint64_t CatchOffset (FILE *file_bai_i, int chr_ID, uint32_t start ){

	int	i,j;
	char	BAMI[4];
	int32_t	n_ref_bai;
	int32_t	n_bin;
	uint32_t	bin;
	int32_t	n_chunk;
//	uint64_t	chunk_beg;
//	uint64_t	chunk_end;
	uint32_t	*n_intv;
	uint64_t	**ioffset;
	
	uint64_t	offset_beg;
//	uint64_t	offset_bgzf;
//	uint64_t	offset_decomp;
	uint32_t	start_index;

	//BAI File
	
	rewind(file_bai_i);
	fread(BAMI,4, sizeof(char), file_bai_i);
	fread(&n_ref_bai, 1, sizeof(int32_t), file_bai_i);
	ioffset	= calloc(sizeof(uint64_t *), n_ref_bai);
	n_intv = calloc(sizeof(int32_t), n_ref_bai);

	for(i = 0;i < n_ref_bai;i++){

		fread(&n_bin,1,sizeof(int32_t),file_bai_i);
		for (j = 0;j < n_bin;j++){
			fread(&bin,1,sizeof(uint32_t),file_bai_i);
			fread(&n_chunk,1,sizeof(int32_t),file_bai_i);
/*
			for(m = 0;m < n_chunk;m++){
				fread(&chunk_beg,1,sizeof(uint64_t),file_bai_i);
				fread(&chunk_end,1,sizeof(uint64_t),file_bai_i);
			}
*/
			fseek(file_bai_i, 2*n_chunk*sizeof(uint64_t), SEEK_CUR);
		}
		fread(&n_intv[i], 1, sizeof(int32_t), file_bai_i);
		printf("%d\t%d\n", i, n_intv[i]);
		if (chr_ID == i){
			ioffset[i] = calloc(sizeof(uint64_t), n_intv[i]);
			fread(ioffset[i], n_intv[i], sizeof(uint64_t), file_bai_i);
		}else {
			fseek(file_bai_i, n_intv[i]*sizeof(uint64_t), SEEK_CUR);
		}	
//		for (j = 0; j <  n_intv[i]; j++){
//			printf("%d\t%lld\t%lld\t%lld\n", j, ioffset[i][j], ioffset[i][j]>>16, ioffset[i][j]&65535);
//		}

		if (chr_ID == i){
		//	printf("%u\t%u\n",ioffset[0], start/16384);
			start_index	= start >> 14;
			if (start_index >= n_intv[i]){
				return 0;			
			}
			offset_beg = ioffset[i][start_index];
			int j;
/*
			offset_bgzf = offset_beg >> 16;
			offset_decomp = offset_beg&65535;
			printf("OFF:%d\t%u\t%lu\t%lu\t%u\t%u\n\n", chr_ID, start, offset_beg, offset_bgzf, offset_decomp, n_intv[i]);
			k = start_index;
			k = (k > 20)?k-20:0;

			for ( ; k <= start/16384; k++){
				offset_bgzf = ioffset[i][k] >> 16;
				offset_decomp = ioffset[i][k]&65535;
				printf("OFF:%d\t%u\t%lu\t%lu\t%u\t%u\n", chr_ID, start, ioffset[i][k], offset_bgzf, offset_decomp, n_intv[i]);
			}
*/
			while ( offset_beg == 0 && (chr_ID != 0 || start_index > 0)){
				if (start_index == 0){
					i--;
					printf("Warning!!\n");
					exit(0);
					if (i < 0 || n_intv[i] == 0){
				//		printf("OFF:%d\t%u\t%lu\t%lu\t%u\t%u\n\n", chr_ID, start, offset_beg, offset_bgzf, offset_decomp, n_intv[i]);
						return 0;
					}
					start_index = n_intv[i] - 1;
				}else {
					start_index --;
				}
				//	printf("OFF:%d\t%u\t%lu\t%lu\t%u\t%u\n\n", chr_ID, start, offset_beg, offset_bgzf, offset_decomp, n_intv[i]);
				offset_beg = ioffset[i][start_index];
			}
			

//			free(ioffset);
			//return offset_decomp;
			for (j = 0; j <=i; j++){
				free(ioffset[j]);
			}
			free(ioffset);
			free(n_intv);

//			offset_bgzf = offset_beg >> 16;
//			offset_decomp = offset_beg&65535;
//			printf("\nOFF:%d\t%u\t%lu\t%lu\t%u\t%u\n\n", chr_ID, start, offset_beg, offset_bgzf, offset_decomp, n_intv[i]);
			return offset_beg;
		}
		if (i == chr_ID){
			break;	
		}
	}
	free(ioffset);
	free(n_intv);
	//Offset Bam File
	return 0;
}


uint8_t *CatchBamHeader (FILE *file_bam_i, bamHeader *BamHeader, uint8_t *stream_i, uint8_t *stream_o, uint8_t *buffer, uint8_t *top){
	int	len_data;
	int	i;	
	int counter;
	int32_t	data_length;
//	int32_t	l_text;
	int32_t l_name;
	uint8_t	*address;
	
	z_stream	infstream;
	bgfzTail	BgfzTail;

	
	len_data = BGFZBlock(file_bam_i);
	decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
	fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);

	memcpy(buffer, stream_o, BgfzTail.I_size);
	address = buffer;	
	
	memcpy(BamHeader->magic, address, sizeof(char)*4);

	memcpy(&(BamHeader->l_text), address+4, sizeof(int32_t));

	BamHeader->text = calloc(BamHeader->l_text + 1,sizeof (char));
//	printf("l_text: %d\n",(BamHeader->l_text));

	address += 8;
	data_length = BgfzTail.I_size - 8;
	while((BamHeader->l_text) - data_length > 0){
		len_data = BGFZBlock(file_bam_i);
		decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
		fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
		memcpy(buffer,stream_o, BgfzTail.I_size);
		address = buffer;	
		
		(BamHeader->l_text) = (BamHeader->l_text) - data_length;
		data_length = BgfzTail.I_size;
	}

	//memcpy(BamHeader->text, address, sizeof(char)*(BamHeader->l_text));
	memcpy(&(BamHeader->n_ref), address + (BamHeader->l_text),sizeof(int32_t));
	address += (4+(BamHeader->l_text));

	BamHeader->chr_name	= malloc(BamHeader->n_ref*sizeof(char *));	
	BamHeader->chr_length	= calloc(BamHeader->n_ref,sizeof(int));

	top = buffer + BgfzTail.I_size;
	counter	= top - address;

	for (i = 0;i < BamHeader->n_ref;i++){
		if (top - address < 4){
			counter = top - address;
			len_data = BGFZBlock(file_bam_i);
			decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
			fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
			//memcpy(buffer, address, sizeof(uint8_t)*counter);
			memmove(buffer, address, sizeof(uint8_t)*counter);

			memcpy(&buffer[counter], stream_o, BgfzTail.I_size);
			top	= buffer + counter + BgfzTail.I_size;
			address = buffer;
		}
		
		memcpy(&l_name, address, sizeof(int32_t));
//		address += 4;
//		chr_name[i]	= calloc(l_name, sizeof(char));
		if (top - address < l_name + 4){
			counter = top - address;
			len_data = BGFZBlock(file_bam_i);
			decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
			fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
			//memcpy(buffer,address,sizeof(uint8_t)*counter);
			memmove(buffer,address,sizeof(uint8_t)*counter);
			memcpy(buffer+counter, stream_o, BgfzTail.I_size);
			top	= buffer + counter + BgfzTail.I_size;
			address = buffer;
		}

//		memcpy(chr_name[i], address, sizeof(char)*l_name);
//		memcpy(&chr_length[i], address+l_name, sizeof(int32_t));
//		address += (l_name+4);
		//printf("%s\t%d\n",chr_name[i], chr_length[i]);
	
//		printf("%d\n", name_len);	
//		BamHeader->chr_name[i]	= calloc(name_len,sizeof(char));
		BamHeader->chr_name[i]	= calloc(l_name,sizeof(char));
		address += refInformation(address, BamHeader->chr_name[i], &(BamHeader->chr_length[i]));
		//printf("%s\t%d\n",BamHeader->chr_name[i], BamHeader->chr_length[i]);
	}
	
	counter	= top - address;
	//counter = BgfzTail.I_size - (address - buffer);
	
	//memcpy(buffer,address,sizeof(uint8_t)*counter);
	memmove(buffer,address,sizeof(uint8_t)*counter);

	top	= buffer + counter;
	//printf("%u\t%u\t%u\n",buffer, address, top);
	
	return top;
}

int SelectReads (uint8_t *stream, alignmentHeader *AlignmentHeader, int QualCut, int NumQualCut){
	
	int	i;
//	char	*read_name;
//	uint32_t *cigar;
//	uint8_t	*seq;
	char	*qual;

	int 	NumQualLess = 0;

//	MemoryCopy(&stream,(void **)&read_name	, AlignmentHeader->l_read_name	, sizeof(char));
//	MemoryCopy(&stream,(void **)&cigar	, AlignmentHeader->n_cigar_op	, sizeof(uint32_t));
//	MemoryCopy(&stream,(void **)&seq	, (AlignmentHeader->l_seq+1)/2	, sizeof(uint8_t));
	stream	+= (AlignmentHeader->l_read_name + (AlignmentHeader->n_cigar_op) * sizeof(uint32_t) +  (AlignmentHeader->l_seq+1)/2);

//	MemoryCopy(&stream,(void **)&qual	, AlignmentHeader->l_seq	, sizeof(char));
	qual	= (char *)stream;	
	
	for (i = 0;i < AlignmentHeader->l_seq; i++){
		if (qual[i] < QualCut){
			NumQualLess++;
		}
	}

//	free(read_name);
//	free(cigar);
//	free(seq);
//	free(qual);

	if (NumQualLess < NumQualCut){
		return 1;
	}else {
		return 0;
	}
}



int reg2bin(int beg, int end){
	--end;
	if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
	if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
	if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
	if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
	if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
	return 0;
}
int	decodeBedFile(FILE *file_bed_i, bedTable *BedTable, bamHeader *BamHeader){

	char	line[LINE_LEN];
	int	num_target	= 0;
	char	chromosome_name[1000];
	int	i;
	int	j	= 0;

	//Bed File
	while	(fgets(line, LINE_LEN, file_bed_i) != NULL){
		num_target++;
	}
	//printf("%d\n",num_target);
	BedTable->table_start	= calloc (BamHeader->n_ref, sizeof(uint32_t));
	BedTable->table_end	= calloc (BamHeader->n_ref, sizeof(uint32_t));
	BedTable->table_max	= calloc (BamHeader->n_ref, sizeof(uint32_t));
	
	BedTable->start	= calloc (num_target, sizeof(uint32_t));
	BedTable->end	= calloc (num_target, sizeof(uint32_t));
	rewind(file_bed_i);

	for(i = 0;i < num_target;i++){
		fscanf(file_bed_i, "%s %d %d", chromosome_name, &BedTable->start[i], &BedTable->end[i]);
		fgets(line, LINE_LEN, file_bed_i);	
	
		if (i == 0){
			for(j = 0; j < BamHeader->n_ref; j++){
				if (strcmp(chromosome_name, BamHeader->chr_name[j])==0){
					BedTable->table_start[j] = i;
					BedTable->table_max[j] = BedTable->end[i];
					break;
				}
			}
		}else if (strcmp(chromosome_name, BamHeader->chr_name[j])!=0){
			BedTable->table_end[j] = i;
			if (BedTable->end[i] > BedTable->table_max[j]){
				BedTable->table_max[j] = BedTable->end[i];
			}
			for(j = 0; j < BamHeader->n_ref; j++){
				if (strcmp(chromosome_name, BamHeader->chr_name[j])==0){
					BedTable->table_start[j] = i;
					BedTable->table_max[j] = BedTable->end[i];
					break;
				}
			}
		}else {
			if (BedTable->end[i] > BedTable->table_max[j]){
				BedTable->table_max[j] = BedTable->end[i];
			}	
		}
		if (j == BamHeader->n_ref){
			printf("Error: %d is over the reference\n", j);
			return -1;
		}
	//printf("%s %d %d\n", chromosome_name, BedTable.start[i], BedTable.end[i]);
	}
	BedTable->table_end[j] = i;
	return num_target;
}

int	decodeBedFile_Merge(FILE *file_bed_i, bedTable *BedTable, bamHeader *BamHeader){

	char	line[LINE_LEN];
	int	num_target	= 0;
	char	chromosome_name[1000];
	int	i;
	int	j = 0;
	int	k;

	//Bed File
	while	(fgets(line, LINE_LEN, file_bed_i) != NULL){
		num_target++;
	}
	//printf("%d\n",num_target);
	BedTable->table_start	= calloc (BamHeader->n_ref, sizeof(uint32_t));
	BedTable->table_end	= calloc (BamHeader->n_ref, sizeof(uint32_t));
	BedTable->table_max	= calloc (BamHeader->n_ref, sizeof(uint32_t));
	
	BedTable->start	= calloc (num_target, sizeof(uint32_t));
	BedTable->end	= calloc (num_target, sizeof(uint32_t));
	rewind(file_bed_i);

	for(i = 0, k = 0;k < num_target;i++, k++){
		fscanf(file_bed_i, "%s %d %d", chromosome_name, &BedTable->start[i], &BedTable->end[i]);
		fgets(line, LINE_LEN, file_bed_i);	
	
		if (i == 0){
			for(j = 0; j < BamHeader->n_ref; j++){
				if (strcmp(chromosome_name, BamHeader->chr_name[j])==0){
					BedTable->table_start[j] = i;
					BedTable->table_max[j] = BedTable->end[i];
					break;
				}
			}
		}else if (strcmp(chromosome_name, BamHeader->chr_name[j])!=0){
			BedTable->table_end[j] = i;
			if (BedTable->end[i] > BedTable->table_max[j]){
				BedTable->table_max[j] = BedTable->end[i];
			}
			for(j = 0; j < BamHeader->n_ref; j++){
				if (strcmp(chromosome_name, BamHeader->chr_name[j])==0){
					BedTable->table_start[j] = i;
					BedTable->table_max[j] = BedTable->end[i];
					break;
				}
			}
		}else {
			if (BedTable->end[i] > BedTable->table_max[j]){
				BedTable->table_max[j] = BedTable->end[i];
			}	
			if (BedTable->start[i] <= BedTable->end[i-1]){
				BedTable->end[i-1] = BedTable->end[i];
				i--;
			}
		}
		if (j == BamHeader->n_ref){
			printf("Error: %d is over the reference\n", j);
			return -1;
		}
	//printf("%s %d %d\n", chromosome_name, BedTable.start[i], BedTable.end[i]);
	}
	BedTable->table_end[j] = i;
	return num_target;
}


int	concernTargetRegion(bedTable *TargetTable, int ref_ID, uint32_t first_pos, uint32_t last_pos){
	int	max_cover_index = -1;
	uint32_t	i;
	int	min_distance	= 999999999;
	int	distance;
	int	start;
	int	end;
	int	min_start	= 999999999;
	int	min_end		= 999999999;

	int	border_margin	= 30;

	for (i =TargetTable->table_start[ref_ID]; i < TargetTable->table_end[ref_ID]; i++){
		if (first_pos < TargetTable->end[i] && last_pos > TargetTable->start[i]){
			if (first_pos >=  TargetTable->start[i]){	distance = first_pos - TargetTable->start[i];}
			else {						distance = TargetTable->start[i] - first_pos;}

			if (last_pos >= TargetTable->end[i]){	distance += (last_pos - TargetTable->end[i]);}
			else {					distance += (TargetTable->end[i] - last_pos);}
			start	= TargetTable->start[i] - first_pos;
			end	= last_pos - TargetTable->end[i];

			if (distance < min_distance || (distance == min_distance && abs(start*end) > abs(min_start*min_end))){
				min_distance = distance;
				min_start	= start;
				min_end	= end;
				max_cover_index = i;
			}
			/*
			if (first_pos >=  TargetTable->start[i]){
				if (last_pos <= TargetTable->end[i]){
					if ( max_cover < 1){
						max_cover = 1;
						max_cover_index = i;
			//			printf("%f\n",1);
					}
				}else {
					if ( max_cover < (float)(TargetTable->end[i]-first_pos) /(last_pos - first_pos)){
						max_cover = (float)(TargetTable->end[i]-first_pos) /(last_pos - first_pos);
						max_cover_index = i;
			//			printf("%f\n",max_cover);
					}
				}
			}else {
				if (last_pos <= TargetTable->end[i]){
					if ( max_cover < (float)(last_pos-TargetTable->start[i])/(last_pos - first_pos)){
						max_cover = (float)(last_pos-TargetTable->start[i])/(last_pos - first_pos);
						max_cover_index = i;
			//			printf("%f\n",max_cover);
					}
				}else {
					if (max_cover < (float)(TargetTable->end[i]-TargetTable->start[i]) /(last_pos - first_pos)){
						max_cover	= (float)(TargetTable->end[i]-TargetTable->start[i]) /(last_pos - first_pos);
						max_cover_index = i;
			//			printf("%f\n",max_cover);
					}
				}
			}
			*/
		}
	}
//	printf("#%u\t%u\t%u\n",min_distance, first_pos, last_pos);
	if (min_start < -border_margin || min_start > border_margin || min_end < -border_margin || min_end > border_margin){
		return -1;	
	}else {
		return max_cover_index;
	}
}
 
void CalculateStatInfo (FILE *file_region_o, bedTable *BedTable, regionInformationLite *RegionInformation, int idx, bamHeader *BamHeader, int ref_ID){
	
	int	j;

	int	length;
	int base;
	int	valid_length;
	double valid_percent;
	float	ave_cov;
	uint32_t	cov_Q2;
	uint32_t	cov_Q3;

	length	= BedTable->end[idx]-BedTable->start[idx];
	base	= 0;
	valid_length	= 0;
	valid_percent	= 0.0;
	ave_cov	= 0.0;
	cov_Q3	= 0;
	cov_Q2	= 0;

	if (RegionInformation[idx].count != 0){
		for (j = 0; j < length; j++){
			if (RegionInformation[idx].ALL[j] > 0){
				base += RegionInformation[idx].ALL[j];
				valid_length++;
			}
		}

		valid_percent	= ((double)valid_length)/length;
		ave_cov	=  ((float)base)/valid_length;
		qsort ( RegionInformation[idx].ALL, length, sizeof(uint32_t), compare_uint32_t);
		cov_Q3	= RegionInformation[idx].ALL[(int)(0.5 + (float)length * 0.75)];
		cov_Q2	= RegionInformation[idx].ALL[(int)(0.5 + (float)length * 0.5)];		
	}

	fprintf(file_region_o, "%s\t", BamHeader->chr_name[ref_ID]);	//chromosome
	fprintf(file_region_o, "%u\t", BedTable->start[idx]);	//start
	fprintf(file_region_o, "%u\t", BedTable->end[idx]);	//end
	fprintf(file_region_o, "%u\t", length);	//length
	fprintf(file_region_o, "%u\t", RegionInformation[idx].count);	//count
	fprintf(file_region_o, "%u\t", length-valid_length);	//length
	fprintf(file_region_o, "%u\t", base);	//length
	fprintf(file_region_o, "%f\t", ave_cov);	//length
	fprintf(file_region_o, "%f\t", valid_percent);	//length
	fprintf(file_region_o, "%u\t", cov_Q2);	//length
	fprintf(file_region_o, "%u\n", cov_Q3);	//length
		
}

uint8_t *CatchBamHeader_v1 (FILE *file_bam_i, bamHeader *BamHeader, uint8_t *stream_i, uint8_t *stream_o){
	int	len_data;
	int	i;	
	int counter;
	int32_t	data_length;
	int32_t l_name;
	uint8_t	*address;
	uint8_t *top;

	z_stream	infstream;
	bgfzTail	BgfzTail;

	
	len_data = BGFZBlock(file_bam_i);
	decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
	fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);

	address = stream_o;	
	
	memcpy(BamHeader->magic, address, sizeof(char)*4);
	memcpy(&(BamHeader->l_text), address+4, sizeof(int32_t));

/*	will not be used temporarily */
//	BamHeader->text = calloc(BamHeader->l_text + 1,sizeof (char));

	address += 8;
	data_length = BgfzTail.I_size - 8;
	while((BamHeader->l_text) - data_length > 0){
		len_data = BGFZBlock(file_bam_i);
		decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
		fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
		address = stream_o;	
		
		(BamHeader->l_text) = (BamHeader->l_text) - data_length;
		data_length = BgfzTail.I_size;
	}

	memcpy(&(BamHeader->n_ref), address + (BamHeader->l_text),sizeof(int32_t));
	address += (4+(BamHeader->l_text));

	BamHeader->chr_name	= malloc(BamHeader->n_ref*sizeof(char *));	
	BamHeader->chr_length	= calloc(BamHeader->n_ref,sizeof(int));

	top	= stream_o + BgfzTail.I_size;
	counter	= top - address;

	for (i = 0;i < BamHeader->n_ref;i++){
		if (top - address < 4){
			counter = top - address;
			len_data = BGFZBlock(file_bam_i);
			memmove(stream_o, address, sizeof(uint8_t)*counter);
			decompressBlock(&infstream, stream_i, stream_o+counter, len_data, file_bam_i);
			fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
			top	= stream_o + counter + BgfzTail.I_size;
			address = stream_o;
		}
		
		memcpy(&l_name, address, sizeof(int32_t));

		if (top - address < l_name + 4){
			counter = top - address;
			len_data = BGFZBlock(file_bam_i);
			memmove(stream_o, address, sizeof(uint8_t)*counter);
			decompressBlock(&infstream, stream_i, stream_o+counter, len_data, file_bam_i);
			fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
			top	= stream_o + counter + BgfzTail.I_size;
			address = stream_o;
		}

	
		BamHeader->chr_name[i]	= calloc(l_name,sizeof(char));
		address += refInformation(address, BamHeader->chr_name[i], &(BamHeader->chr_length[i]));
		//printf("%s\t%d\n",BamHeader->chr_name[i], BamHeader->chr_length[i]);
	}
	
	counter	= top - address;
	memmove( stream_o, address, sizeof(uint8_t)*counter);
	top	= stream_o + counter;
	
	return top;
}

#endif
