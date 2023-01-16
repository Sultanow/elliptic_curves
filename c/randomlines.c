#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main() {
    FILE *fp;
    char command[100];
    char line[1000];
    int lines;
    int random_lines[100];
    int i;
    
    // Get the number of lines in the file using the wc -l command
    sprintf(command, "wc -l textfile.txt");
    fp = popen(command, "r");
    fscanf(fp, "%d", &lines);
    pclose(fp);

    // Seed the random number generator
    srand(time(0));

    // Generate 100 random line numbers
    for (i = 0; i < 100; i++) {
        random_lines[i] = rand() % lines + 1;
    }

    // Print out the 100 random lines from the file
    for (i = 0; i < 100; i++) {
        sprintf(command, "sed -n '%dp' textfile.txt", random_lines[i]);
        fp = popen(command, "r");
        while (fgets(line, sizeof(line), fp) != NULL) {
            printf("%s", line);
        }
        pclose(fp);
    }

    return 0;
}
