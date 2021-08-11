#ifndef FILE_IO_H
#define FILE_IO_H

#include <stdio.h>
#include <stdlib.h>

#define OVERWRITE 1
#define APPEND 2

void file_writer(char fname[], char msg[], int WRITE_TYPE);
char *file_reader(char fname[]);


int test_file_io_functionality();

#endif