#include <stdlib.h>
#include <stdio.h>
#include "htslib/vcf.h"
#include <assert.h>
#include <stdint.h>
#include <sys/time.h>

int read_multiple_lines(htsFile* fptr, bcf_hdr_t* hdr, bcf1_t* array[], unsigned max)
{
    int i = 0;
    int status = 0;
    for(i=0;i<max;++i)
    {
        status = bcf_read(fptr, hdr, array[i]);
        if(status < 0)
            break;
    }
    return i;
}

#define SIZE 10000
void speed_check(char** argv)
{
    int i = 0;
    int j = 0;
    //open file and read record
    htsFile* fptr = hts_open(argv[1], "r");
    assert(fptr);
    bcf_hdr_t* hdr = bcf_hdr_read(fptr);
    assert(hdr);
    
    int info_ids[50];
    int num_info_fields = 0;
    for(i=0;i<hdr->n[BCF_DT_ID];++i)
        if(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, i))
            info_ids[num_info_fields++] = i;
    printf("Num info fields %d\n",num_info_fields);
    for(i=0;i<1000;++i)
    {
        int x = rand()%num_info_fields;
        int y = rand()%num_info_fields;
        int tmp = info_ids[x];
        info_ids[x] = info_ids[y];
        info_ids[y] = tmp;
    }
    bcf1_t* array[SIZE];
    for(i=0;i<SIZE;++i)
    {
      array[i] = bcf_init();
    }
    uint64_t sum_num_info = 0;
    unsigned num_lines = 0;
    struct timeval t1, t2;
    uint64_t total_time = 0;
    while(1)
    {
        int num = read_multiple_lines(fptr, hdr, array, SIZE);
        if(num == 0)
            break;
        for(i=0;i<num;++i)
            bcf_unpack(array[i], BCF_UN_INFO);
        gettimeofday(&t1, 0);
        {
            for(i=0;i<num;++i)
            {
                if(array[i]->n_info > 1)
                {
                    int k;
#define REPEAT_TIMES 1000
                    for(k=0;k<REPEAT_TIMES;++k)
                    {
                        for(j=0;j<num_info_fields;++j)
                            sum_num_info += (bcf_get_info_idx(array[i], info_ids[j])) ? 1 : 0;
                        ++num_lines;
                    }
                }
            }
        }
        gettimeofday(&t2, 0);
        total_time += ((t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec));
    }
    fprintf(stderr,"Num info fields %lu total num lines %d repeats %d time %lf\n",sum_num_info, num_lines, REPEAT_TIMES, ((double)total_time)*1e-6);
    for(i=0;i<SIZE;++i)
        bcf_destroy(array[i]);
    bcf_hdr_destroy(hdr);
    hts_close(fptr);
}

void functional_check(char** argv)
{
    //open file and read record
    htsFile* fptr = hts_open(argv[1], "r");
    assert(fptr);
    bcf_hdr_t* hdr = bcf_hdr_read(fptr);
    assert(hdr);
    bcf1_t* line = 0;

    int status = 0;
    line = bcf_init();
    int i = 0;
    int info_ids[50];
    int num_info_fields = 0;
    for(i=0;i<hdr->n[BCF_DT_ID];++i)
        if(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, i))
            info_ids[num_info_fields++] = i;
#define ERASE_SIZE 2
#define ADD_SIZE 1
    int erase_ids[ERASE_SIZE];
    for(i=0;i<ERASE_SIZE;++i)
    {
        erase_ids[i] = info_ids[rand()%num_info_fields];
        fprintf(stderr,"Will erase from all bcf1_t : %s\n",bcf_hdr_int2id(hdr, BCF_DT_ID, erase_ids[i]));
    }
    int add_ids[ADD_SIZE];
    for(i=0;i<ADD_SIZE;++i)
    {
        add_ids[i] = info_ids[rand()%num_info_fields];
        fprintf(stderr,"Will add to all bcf1_t : %s\n",bcf_hdr_int2id(hdr, BCF_DT_ID, add_ids[i]));
    }

    unsigned data_array[] = { 56, 56, 56, 56, 67, 64, 60 };

    /*int hdr_length = 0;*/
    /*printf("%s",bcf_hdr_fmt_text(hdr, 0, &hdr_length));*/

    unsigned num_lines = 0;
    kstring_t debug_string = { 0, 0, 0 };
    int nbuffer = 1024;
    int* buffer = (int*)malloc(nbuffer*sizeof(int));
    assert(buffer != 0);
    while(1)
    {
        status = bcf_read(fptr, hdr, line);
        if(status < 0)
            break;
        bcf_unpack(line, BCF_UN_INFO);

        for(i=0;i<ERASE_SIZE;++i)
            bcf_update_info(hdr, line, bcf_hdr_int2id(hdr, BCF_DT_ID, erase_ids[i]), 0, 0, bcf_hdr_id2type(hdr, BCF_HL_INFO, erase_ids[i])); 
        for(i=0;i<ADD_SIZE;++i)
            bcf_update_info(hdr, line, bcf_hdr_int2id(hdr, BCF_DT_ID, add_ids[i]), (void*)data_array, 
                   bcf_hdr_id2number(hdr, BCF_HL_INFO, add_ids[i]), bcf_hdr_id2type(hdr, BCF_HL_INFO, add_ids[i])); 

	for(i=0;i<ERASE_SIZE;++i)
	{
	  const char* key = bcf_hdr_int2id(hdr, BCF_DT_ID, erase_ids[i]);
	  assert(bcf_get_info_values(hdr, line, key, (void**)&buffer, &nbuffer, bcf_hdr_id2type(hdr, BCF_HL_INFO, erase_ids[i])) == -3
	      && "Found erased field in bcf1_t structure");
	}
	for(i=0;i<ADD_SIZE;++i)
	{
	  const char* key = bcf_hdr_int2id(hdr, BCF_DT_ID, add_ids[i]);
	  int ht_type = bcf_hdr_id2type(hdr, BCF_HL_INFO, add_ids[i]);
	  int num_values = bcf_get_info_values(hdr, line, key, (void**)&buffer, &nbuffer, ht_type);
	  assert(num_values >= 0 && "Added field not found in bcf1_t structure");
	  if(ht_type != BCF_HT_FLAG)
	  {
	    unsigned element_size = bcf_hdr_id2type(hdr, BCF_HL_INFO, add_ids[i]);
	    unsigned j = 0u;
	    for(j=0;j<num_values*element_size;++j)
	      assert(((char*)buffer)[j] == ((char*)data_array)[j] && "Added field value differs");
	  }
	}

        /*bcf_update_info(hdr, line, "END", &new_value, 1, BCF_HT_INT);*/
        //calls bcf1_sync internally
        debug_string.l = 0;
	/*vcf_format(hdr, line, &debug_string);*/
	/*printf("%s",debug_string.s);*/
        ++num_lines;
    }
    free(buffer);
    if(debug_string.m > 0)
        free(debug_string.s);
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    hts_close(fptr);
#undef ERASE_SIZE
#undef ADD_SIZE
}

void translate_check(char** argv)
{
    //open file and read record
    htsFile* fptr = hts_open(argv[1], "r");
    assert(fptr);
    bcf_hdr_t* hdr = bcf_hdr_read(fptr);
    assert(hdr);
    //duplicate
    bcf_hdr_t* dup_hdr = bcf_hdr_dup(hdr);
    assert(dup_hdr);

    bcf1_t* line = 0;

    int status = 0;
    line = bcf_init();
    int i = 0;
    int info_ids[50];
    int num_info_fields = 0;
    for(i=0;i<hdr->n[BCF_DT_ID];++i)
        if(bcf_hdr_idinfo_exists(hdr, BCF_HL_INFO, i))
            info_ids[num_info_fields++] = i;
#define ERASE_SIZE 1
#define ADD_SIZE 1
    int erase_ids[ERASE_SIZE];
    for(i=0;i<ERASE_SIZE;++i)
    {
        erase_ids[i] = info_ids[rand()%num_info_fields];
        fprintf(stderr,"Removing from header : %s\n",bcf_hdr_int2id(hdr, BCF_DT_ID, erase_ids[i]));
    }
    int add_ids[ADD_SIZE];
    for(i=0;i<ADD_SIZE;++i)
    {
        add_ids[i] = info_ids[rand()%num_info_fields];
        fprintf(stderr,"Re-adding to header : %s\n",bcf_hdr_int2id(hdr, BCF_DT_ID, add_ids[i]));
    }
    for(i=0;i<ERASE_SIZE;++i)
    {
        bcf_hdr_remove(dup_hdr, BCF_HL_INFO, bcf_hdr_int2id(hdr, BCF_DT_ID, erase_ids[i]));
        /*bcf_hrec_t* del_hrec = bcf_hdr_get_hrec(dup_hdr, BCF_HL_INFO, "ID", bcf_hdr_int2id(hdr, BCF_DT_ID, erase_ids[i]), 0);*/
        /*assert(del_hrec);*/
    }
    for(i=0;i<ADD_SIZE;++i)
    {
        bcf_hrec_t* del_hrec = bcf_hdr_get_hrec(dup_hdr, BCF_HL_INFO, "ID", bcf_hdr_int2id(hdr, BCF_DT_ID, add_ids[i]), 0);
        assert(del_hrec);
        bcf_hrec_t* readd = bcf_hrec_dup(del_hrec);
        bcf_hdr_remove(dup_hdr, BCF_HL_INFO, bcf_hdr_int2id(hdr, BCF_DT_ID, add_ids[i]));
        bcf_hdr_add_hrec(dup_hdr, readd);
    }

    /*int hdr_length = 0;*/
    /*printf("%s",bcf_hdr_fmt_text(dup_hdr, 0, &hdr_length));*/

    unsigned num_lines = 0;
    kstring_t debug_string = { 0, 0, 0 };
    int nbuffer = 1024;
    int* buffer = (int*)malloc(nbuffer*sizeof(int));
    assert(buffer != 0);
    while(1)
    {
        status = bcf_read(fptr, hdr, line);
        if(status < 0)
            break;
        bcf_unpack(line, BCF_UN_INFO);
        bcf_translate(dup_hdr, hdr, line);
	
	for(i=0;i<ERASE_SIZE;++i)
	{
	  const char* key = bcf_hdr_int2id(hdr, BCF_DT_ID, erase_ids[i]);
	  assert(bcf_get_info_values(hdr, line, key, (void**)&buffer, &nbuffer, bcf_hdr_id2type(hdr, BCF_HL_INFO, erase_ids[i])) == -3
	      && "Found erased field in bcf1_t structure");
	}

        debug_string.l = 0;
	/*vcf_format(dup_hdr, line, &debug_string);*/
	/*printf("%s",debug_string.s);*/
        ++num_lines;
    }
    free(buffer);

    if(debug_string.m > 0)
        free(debug_string.s);
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    hts_close(fptr);
}

int main(int argc, char** argv)
{
  assert(argc >= 2 && "Requires 1 argument <vcf/bcf> file name");
  fprintf(stderr,"Running functional check for direct access ..\n");
  functional_check(argv);
  fprintf(stderr,"Running bcf_translate check for direct access ..\n");
  translate_check(argv);
  fprintf(stderr,"Running speed check for direct access ..\n");
  speed_check(argv);

  return 0;
}

