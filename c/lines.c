#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    FILE *file = fopen("candidates.txt", "r");
    if (file == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }

    int data[5];
    char string_example[1000];
    while (fscanf(file, "[ %*d, %d, %d, %d, %d, %d, %*d ]\n", &data[0], &data[1], &data[2], &data[3], &data[4]) == 5) {
    
     sprintf(string_example, "%d, %d, %d, %d, %d", data[0], data[1], data[2], data[3], data[4]);
    
     char command[1000];
     snprintf(command, 1000, "./lpdata0 test [%s] 2e17", string_example);
     printf("%s\n",command);
      // Make the system call
     system(command);
       
    }

    fclose(file);
    return 0;
}
