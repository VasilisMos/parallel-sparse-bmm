#include "file_io.hpp"

// TODO
char *file_reader(char fname[]){
    FILE *f;

    f = fopen(fname, "a");

    if (f == NULL){
        printf("Couldn't Open %s, Exiting\n",fname);
        exit(1);
    }

    fclose(f);
}

void file_writer(char fname[], char msg[], int WRITE_TYPE){
    FILE *f;

    if(WRITE_TYPE == APPEND)
        f = fopen(fname, "a");
    
    if(WRITE_TYPE == OVERWRITE)
        f = fopen(fname, "w");

    if (f == NULL){
        printf("Couldn't Open %s, Exiting\n",fname);
        exit(1);
    }
        
    fprintf(f, "%s\n", msg);

    fclose(f);
}

int test_file_io_functionality(){
    char test_msg[] = "Hello World";
    char fname[] = "test_write_file.txt";

    file_writer(fname, test_msg, APPEND);


    return 0;

}