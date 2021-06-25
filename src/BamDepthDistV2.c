#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include "zlib.h"

#include "BamCommonLibrary.h"

int	BamDepthDist(FILE *file_bam_i, FILE *file_bai_i, FILE *file_bed_i, toolsFlags *ToolsFlags){

	uint8_t *stream_i;//[65536];
	uint8_t *stream_o;//[65536];
	uint8_t	*top;

	int	j;
	uint32_t n, k;
	int	len_data;
	int	counter;
	int	ref_ID	= -1;

	z_stream infstream;
	bgfzTail BgfzTail;
	
	alignmentHeader	AlignmentHeader;

	uint8_t *address;

	uint64_t	offset_beg;
	uint64_t	offset_bgzf;
	uint64_t	_offset_bgzf;
	uint64_t	offset_decomp;
	posCoverage 	*PosCoverage	= NULL;
	bamHeader	BamHeader;
	bedTable	BedTable;
	bedTable	TargetTable;
	baiTable	BaiTable;
	
	int32_t	end_position	= 0;


	uint32_t	max_position;
	uint64_t	Total_Depth = 0;
	uint64_t	temp_Depth	= 0;
	uint64_t	Total_Length = 0;
	uint64_t	Cover_Length = 0;
	int32_t	max_cover_index;

	uint64_t	*threshold	= NULL;
	uint64_t	*num_coverage_chr	= NULL;
	uint64_t	*num_coverage_all	= NULL;
	uint64_t	coverage;
	
	int	num_threshold = 0;
	int	index_threshold = 0;
	char	*token;
	char	*duplicate_string;

	uint32_t	region_start	= 0;
	uint32_t	region_end	= 0;
	uint32_t	region_offset	= 0;
	
	uint32_t	end_read_pos	= 0;
	uint32_t	max_base_pos	= 0;


	if (ToolsFlags->flag_coverage == 1){
		duplicate_string = strdup(ToolsFlags->coverage_threshold);
		//token = strtok( ToolsFlags->coverage_threshold, ",");
		token = strtok( duplicate_string, ",");
		while( token != NULL ){
			num_threshold++;
			/* Get next token: */
			token = strtok( NULL, ",");
		}

		threshold = calloc(num_threshold, sizeof(int));
		num_coverage_all = calloc(num_threshold, sizeof(uint64_t));
		num_coverage_chr = calloc(num_threshold, sizeof(uint64_t));

		num_threshold	= 0;
		//printf("%s\n", ToolsFlags->coverage_threshold);
		token = strtok( ToolsFlags->coverage_threshold, ",");
		while( token != NULL ){
			/* While there are tokens in "string" */
			threshold[num_threshold] = atoi(token);
//			printf("%d\n", atoi(token));
			token = strtok( NULL, ",");
			if (threshold[num_threshold] == 0){
				printf("Error: value of threshold\n");
				return -1;
			}
			num_threshold++;
			/* Get next token: */
		}
//		printf("%d\n",num_threshold);
		printf("THRESH");
		for ( index_threshold = 0; index_threshold < num_threshold; index_threshold++){
			printf("\t%luX", threshold[index_threshold]);
		}
		printf("\t** (THRESH: threshold)\n");
		
	}

	

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
	if (ToolsFlags->flag_hide == 1){
		for (i = 0;i < BamHeader.n_ref;i++){
			printf("%s\t%d\n",BamHeader.chr_name[i], BamHeader.chr_length[i]);
		}
	}
	*/

	if (ToolsFlags->flag_target){
		decodeBedFile(ToolsFlags->file_target, &TargetTable, &BamHeader);
	}

	createRegionTable(file_bed_i, &BedTable, &BamHeader, ToolsFlags);


	if (ToolsFlags->flag_coverage == 0 && ToolsFlags->flag_simple == 0){
		printf("#CHR\tPOS\tA\tC\tG\tT\tN\tDEL\tTOTAL\n");
	}


	// reference
	//	// bed
	//	//	//	while
		
	BuildBaiTable (file_bai_i, &BaiTable);
	//BAI File
	for (j = 0; j < BamHeader.n_ref; j++){
		ref_ID = j;
		if (BedTable.table_end[ref_ID] == BedTable.table_start[ref_ID]){	// Check whether any region is in this ref_ID. 
			continue;
		}

		max_position = 0;
		end_read_pos	= BamHeader.chr_length[j] + 1;
		max_base_pos	= 0;
		counter = 0;
//		printf("YCL REF:%d\n", j);
		_offset_bgzf	= 0;
		for (k = BedTable.table_start[ref_ID]; k < BedTable.table_end[ref_ID]; k++){
			offset_beg = ReturnOffset ( &BaiTable, ref_ID, BedTable.start[k]);
			//offset_beg = CatchOffset (file_bai_i, ref_ID, BedTable.start[k]);
			offset_bgzf = offset_beg >> 16;
			offset_decomp = offset_beg&65535;

			/* Modified by Yu-Cheng Li 202012 [Start]*/
			// PosCoverage = calloc(BamHeader.chr_length[ref_ID],sizeof(posCoverage));
			region_offset	= BedTable.start[k];
			region_start	= BedTable.start[k];
			region_end	= BedTable.end[k];

			PosCoverage = calloc(BedTable.end[k] - BedTable.start[k],sizeof(posCoverage));
			/* Modified by Yu-Cheng Li 202012 [End]*/
	
//			printf("YCL %lu\t%lu\t%lu\n",offset_beg, offset_bgzf, offset_decomp);			
//
			// Added by Yu-Cheng Li 202102 [Start]
//			printf("YCL ERP\t%d\t%lu\t%lu\n", max_base_pos, offset_bgzf, _offset_bgzf);
			if ((k != BedTable.table_start[ref_ID]) && (max_base_pos < region_start) && (offset_bgzf == _offset_bgzf)){
				offset_beg = 0;
			}
			//

			if (offset_beg == 0){

			}else {
			//Offset Bam File
				_offset_bgzf	= offset_bgzf;
				fseek(file_bam_i,offset_bgzf,SEEK_SET);
				len_data = BGFZBlock(file_bam_i);
				decompressBlock(&infstream, stream_i, stream_o, len_data, file_bam_i);
//				printf("decompress 1 YCL %d\t%d\n", k, offset_bgzf);
				fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);
				top = stream_o + BgfzTail.I_size;
				address = stream_o + offset_decomp;
				counter = BgfzTail.I_size - offset_decomp;
			}
			//Start Unzip and Process Bam File	
			while ( (len_data = BGFZBlock(file_bam_i)) > 0 ){
				if (counter < 65536){
					memmove( stream_o, address, sizeof(uint8_t)*counter);
					decompressBlock(&infstream, stream_i, stream_o+counter, len_data, file_bam_i);
					
					fread(&BgfzTail	,sizeof(bgfzTail),1,file_bam_i);

					top	= stream_o + counter + BgfzTail.I_size;
					address = stream_o;
//					printf("Length YCL %d\t%d\t%d\n", k, BgfzTail.I_size, top-address);

				}else {
					fseek(file_bam_i, -sizeof(bgfzHeader), SEEK_CUR);
				}
				if (top - address > (intptr_t) sizeof(alignmentHeader)){
					memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
				}
				//printf("%u\t%u\t%u\n", top, address, AlignmentHeader.refID);
				while ( top - address >= (AlignmentHeader.block_size + 4) ){
					if (AlignmentHeader.refID > ref_ID || AlignmentHeader.refID < (ref_ID - 1)){
						break;
					}
					//if (AlignmentHeader.refID == ref_ID && AlignmentHeader.pos > BedTable.table_max[ref_ID]){
					if (AlignmentHeader.refID == ref_ID && (uint32_t) AlignmentHeader.pos > BedTable.end[k]){
						break;	
					} 
//					printf("Alignment YCL %d\t%d\t%d\t%d\n", k, AlignmentHeader.block_size, AlignmentHeader.pos, top - address);
					address += sizeof(alignmentHeader);
					//Core_START

					if (AlignmentHeader.refID == ref_ID){
						if ((AlignmentHeader.FLAG&1024)==0 || ToolsFlags->flag_dup==0 ){
							if (!ToolsFlags->flag_filter || SelectReads(address, &AlignmentHeader, 30, 30)){
								if (ToolsFlags->flag_target){
									end_position	= alignment_EndPosition(address, &AlignmentHeader, BamHeader.chr_length[ref_ID]);
//									printf("YCL %u\n", end_position);
									max_cover_index	= concernTargetRegion( &TargetTable, ref_ID, AlignmentHeader.pos, end_position);
									//printf("%u\t%u\t%d\n", AlignmentHeader.pos, end_position, max_cover_index);
									if (max_cover_index >= 0){
										region_start	= BedTable.start[k];
										region_end	= BedTable.end[k];
										if (region_start < TargetTable.start[max_cover_index]){	region_start	=	TargetTable.start[max_cover_index];}
										if (region_end > TargetTable.end[max_cover_index]){	region_end	=	TargetTable.end[max_cover_index];}

										alignment_DepthDist_TR(address, &AlignmentHeader, PosCoverage, region_start, region_end, region_offset);	
									}
								}else {
									end_read_pos = alignment_DepthDist_TR(address, &AlignmentHeader, PosCoverage, region_start, region_end, region_offset);
									if (end_read_pos > max_base_pos){
										max_base_pos	= end_read_pos;
										//printf("YCL MAX %d\t%d\n",AlignmentHeader.pos, max_base_pos);
									}
//									printf("Read alignment YCL %d\n",k);
								}
							}
						}
					}
					//Core_END
			
					address += (AlignmentHeader.block_size + 4) - sizeof(alignmentHeader);

					if (top - address > (intptr_t) sizeof(alignmentHeader)){
						memcpy(&AlignmentHeader, address, sizeof(alignmentHeader));
					}else {
						break;	
					}
				}	// End of 'while' for 'alignment'

				counter = top - address;

			//	printf("%u\t%u\t%u\n", BedTable.start[BedTable.table_start[ref_ID]], AlignmentHeader.pos, BedTable.table_max[ref_ID]);
				if (AlignmentHeader.refID > ref_ID || AlignmentHeader.refID < (ref_ID - 1)){
		//		if (AlignmentHeader.refID > ref_ID){
					break;
				}
				//if (AlignmentHeader.refID == ref_ID && AlignmentHeader.pos > BedTable.table_max[ref_ID]){
				if (AlignmentHeader.refID == ref_ID && (uint32_t) AlignmentHeader.pos > BedTable.end[k]){
					break;
				}
			}	// End of 'while' for 'decompress block'

			if (ToolsFlags->flag_coverage == 1){
				if ( BedTable.end[k] > max_position){
					if (BedTable.start[k] < max_position){
						BedTable.start[k] = max_position;
					}

					for (n = BedTable.start[k]-region_offset; n < BedTable.end[k]-region_offset;n++){
						coverage = (PosCoverage[n].A +
							PosCoverage[n].C +
							PosCoverage[n].G +
							PosCoverage[n].T +
							PosCoverage[n].N +
							PosCoverage[n].DEL);
						Total_Length++;
						for ( index_threshold = 0; index_threshold < num_threshold; index_threshold++){
							if (coverage >= threshold[index_threshold]){	
								num_coverage_chr[index_threshold]++;	
								num_coverage_all[index_threshold]++;	
							}
						}
					}
					max_position = BedTable.end[k];
				}	
			}else if (ToolsFlags->flag_simple == 1){
				if ( BedTable.end[k] > max_position){
					if (BedTable.start[k] < max_position){
						BedTable.start[k] = max_position;
					}
					
					for (n = BedTable.start[k]-region_offset; n < BedTable.end[k]-region_offset;n++){
						temp_Depth = (PosCoverage[n].A +
						PosCoverage[n].C +
						PosCoverage[n].G +
						PosCoverage[n].T +
						PosCoverage[n].N +
						PosCoverage[n].DEL);
						if (temp_Depth > 0){
							Cover_Length++;
							Total_Depth += temp_Depth;
						}
						Total_Length++;
					}
					max_position = BedTable.end[k];	
				}
			}else {
				if ( BedTable.end[k] > max_position){
					if (BedTable.start[k] < max_position){
						BedTable.start[k] = max_position;
					}
					for (n = BedTable.start[k]-region_offset; n < BedTable.end[k]-region_offset;n++){
/*
						printf("%s\t",BamHeader.chr_name[j]);
						printf("%u\t",n + 1 + region_offset);
						printf("%u\t",PosCoverage[n].A);
						printf("%u\t",PosCoverage[n].C);
						printf("%u\t",PosCoverage[n].G);
						printf("%u\t",PosCoverage[n].T);
						printf("%u\t",PosCoverage[n].N);
						printf("%u\t",PosCoverage[n].DEL);
						printf("%u",	PosCoverage[n].A +
							PosCoverage[n].C +
							PosCoverage[n].G +
							PosCoverage[n].T +
							PosCoverage[n].N +
							PosCoverage[n].DEL);					
						printf("\n");
*/
///*
						printf("%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n", 
							BamHeader.chr_name[j],
							n + 1 + region_offset,
							PosCoverage[n].A,
							PosCoverage[n].C,
							PosCoverage[n].G,
							PosCoverage[n].T,
							PosCoverage[n].N,
							PosCoverage[n].DEL,
							PosCoverage[n].A+PosCoverage[n].C+PosCoverage[n].G+PosCoverage[n].T+PosCoverage[n].N+PosCoverage[n].DEL);
//*/
					}
					max_position = BedTable.end[k];	
				}
			}
			free(PosCoverage);

		}
		if (ToolsFlags->flag_hide == 1){
			printf("[Bam File Unzip %d / %d ] %s done\n",ref_ID+1, BamHeader.n_ref, BamHeader.chr_name[ref_ID]);
		}			

		if (ToolsFlags->flag_coverage == 1){
			printf("%s",BamHeader.chr_name[ref_ID]);
			for ( index_threshold = 0; index_threshold < num_threshold; index_threshold++){
				printf("\t%lu",num_coverage_chr[index_threshold]);
				num_coverage_chr[index_threshold] = 0;
			}
			printf("\n");
		}		
	}
	
	if (ToolsFlags->flag_coverage == 1){
		printf("TOT.DP");
		for ( index_threshold = 0; index_threshold < num_threshold; index_threshold++){
			printf("\t%lu",num_coverage_all[index_threshold]);
		}
		printf("\n");
		printf("COV.RT");
		for ( index_threshold = 0; index_threshold < num_threshold; index_threshold++){
			printf("\t%.5lf", (double)num_coverage_all[index_threshold] /Total_Length);
		}
		printf("\n");

	}else if (ToolsFlags->flag_simple == 1){
		printf("#AVE.DP\tCOV.RT\t**(AVE.DP: average depth; COV.RT: coverage rate)\n");
		printf("%.4f\t%.4f\n", (double)Total_Depth / Cover_Length, (double)Cover_Length / Total_Length);
	}
	free(stream_i);	
	free(stream_o);	
	return 0;
}

