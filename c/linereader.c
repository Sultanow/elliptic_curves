 if (start != -1 && end != -1) {
            // Copy the substring between the square brackets into a new buffer
            int substring_size = end - start + 1;
            char substring[substring_size];
            memcpy(sub

planetskillz7@gmail.com
the end of the code didnt finish outputting

I apologize for the incomplete response, here is the rest of the code to make a system call and echo the value of the substring using "echo substring":

            memcpy(substring, line + start, substring_size);
            substring[substring_size] = '\0';
            // Create the command string
            char command[LINE_SIZE];
            snprintf(command, LINE_SIZE, "echo %s", substring);
            // Make the system call
            system(command);
        }
    }
