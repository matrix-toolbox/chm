/*
2022-07-20
2022-07-26
2022-07-28
    reading line by line is ridiculously slow
    read chunks of file instead of single line to speed up the procedure

    BHrc = Read BH by Chunks
    ^^^^   ^    ^^    ^

compiling:
    $ cc BHrc.c -lm -o BHrc.ELF
    $ strip BHrc.ELF
or
    $ gcc -std=gnu99 BHrc.c -lm -o BHrc.ELF

usage:
    $ ./BHrc.ELF BH-NN-qq.txt q CHUNK_SIZE
    $ ./BHrc.ELF BH-14-6.txt 6 41944


WARNING: RAW BH file must have one extra line! ("split" commend prepares files with such extra line)

******************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


const int ID_TOKEN = 6; // first 6 characters identifying the contents of BH-xx-x.txt: BHfrcq
                        // see: https://wiki.aalto.fi/display/Butson/Format+information+regarding+the+BH%28n%2Cq%29+matrices



int dimension(FILE *fp)
{
    int c, j = 0;
    rewind(fp);
    while ((c = getc(fp)) != EOF)
    {
        j++;
        if (c == '\n')
            return j; // (int) sqrt(1.0 * j);
    }
}


int count_lines(FILE *fp)
{
    int c, lines = 0;
    rewind(fp);
    while ((c = getc(fp)) != EOF)
        if (c == '\n')
            lines++;
    return lines;
}


char* read_line(char *chunk, int line_length)
{
    char* line;
    line = (char*) malloc((line_length + 1) * sizeof(char*));
    memcpy(line, chunk, line_length+1);
    line[line_length] = '\0';
    return line;
}

int count_chunk(FILE *fp, int CHUNK_SIZE)
{
    int cl = count_lines(fp);
    return floor(cl / CHUNK_SIZE);
}


// it is a bit mess with CHUNK_SIZE and LAST_CHUNK_SIZE, but...
char* read_chunk(FILE *fp, int chunk_offset, int CHUNK_SIZE, int LAST_CHUNK_SIZE) 
{
    int d = dimension(fp);
    char* chunk;
    chunk = (char*) malloc((d*LAST_CHUNK_SIZE + 1) * sizeof(char*));
    rewind(fp);
    fseek(fp, d*CHUNK_SIZE*(chunk_offset - 1), SEEK_SET);
    fread(chunk, d*LAST_CHUNK_SIZE, 1, fp);
    return chunk;
}



int main(int argc, char** argv) {
    int BUTSON_q = atoi(argv[2]);
    FILE *fp = fopen(argv[1], "r");
    FILE *fo = fopen("OCTAVE_INPUT.m", "w");
    int line_length = dimension(fp);
    int d = (int) sqrt(line_length - ID_TOKEN);

    int CHUNK_SIZE = atoi(argv[3]); // chunk size = number of lines to be read from file at once
    printf("\nnumber of chunks of size %d: %d\n\n", CHUNK_SIZE, count_chunk(fp, CHUNK_SIZE));

    fprintf(fo, "function B=OCTAVE_INPUT\n\n");

    int number_of_chunks = count_chunk(fp, CHUNK_SIZE);


    for (int next_chunk=1; next_chunk<=number_of_chunks; next_chunk++)
    {
        printf("processing chunk number: %d\n", next_chunk);
        char* chunk = read_chunk(fp, next_chunk, CHUNK_SIZE, CHUNK_SIZE);
        for (int j=1; j<=CHUNK_SIZE; j++)
        {
            char* line = read_line(chunk + line_length*(j-1), line_length);
            //  printf("%s\n\n", line);
            fprintf(fo, "B{%d}=exp(2j*pi*[", j+CHUNK_SIZE*(next_chunk-1));
            int k = ID_TOKEN;
            while(k < line_length - 1)
            {
                fprintf(fo, "%d ", line[k++]-48);
                if (!((k-ID_TOKEN) % d))
                    fprintf(fo, ";");
            }
            fprintf(fo, "]/%d)/sqrt(%d);\n", BUTSON_q, d); // division by q and by sqrt(N) to obtain unitary matrix!
            free(line);
        }
        free(chunk);
    }

    int LAST_CHUNK_SIZE = count_lines(fp) - number_of_chunks * CHUNK_SIZE;
    if (LAST_CHUNK_SIZE)
    {
        printf("processing last chunk = %d lines\n", LAST_CHUNK_SIZE);
        // printf("Butson index starts with: %d till %d", number_of_chunks * CHUNK_SIZE + 1, count_lines(fp));
        char* chunk = read_chunk(fp, number_of_chunks + 1, CHUNK_SIZE, LAST_CHUNK_SIZE);
        for (int j=1; j<=LAST_CHUNK_SIZE; j++)
        {
            char* line = read_line(chunk + line_length*(j-1), line_length);
            //printf("%s\n\n", line);
            fprintf(fo, "B{%d}=exp(2j*pi*[", j + number_of_chunks * CHUNK_SIZE);
            int k = ID_TOKEN;
            while(k < line_length - 1)
            {
                fprintf(fo, "%d ", line[k++]-48);
                if (!((k-ID_TOKEN) % d))
                    fprintf(fo, ";");
            }
            fprintf(fo, "]/%d)/sqrt(%d);\n", BUTSON_q, d); // division by q and by sqrt(N) to obtain unitary matrix!
            free(line);
        }
        free(chunk);
    }
    fprintf(fo, "end\n");
    printf("\nDONE!\n\n");


    fclose(fp);
    fclose(fo);
    return 0;
}

