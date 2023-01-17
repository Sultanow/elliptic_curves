#include <stdio.h>
#include <string.h>

#define LINE_SIZE 1000

int main() {
    // Open the text file
    FILE *fp = fopen("file.txt", "r");
    if (fp == NULL) {
        printf("Error opening file\n");
        return 1;
    }

    // Read the file line by line
    char line[LINE_SIZE];
    while (fgets(line, LINE_SIZE, fp) != NULL) {
        // Find the start and end indices of the inner square brackets
        int start = -1, end = -1;
        for (int i = 0; i < strlen(line); i++) {
            if (line[i] == '[') {
                start = i;
            }
            else if (line[i] == ']') {
                end = i;
                break;
            }
        }
        if (start != -1 && end != -1) {
            // Copy the substring between the square brackets into a new buffer
            int substring_size = end - start + 1;
            char substring[substring_size];
            memcpy(substring, line + start, substring_size);
            substring[substring_size] = '\0';
            printf("%s\n", substring);
        }
    }

    // Close the file
    fclose(fp);

    return 0;
}
