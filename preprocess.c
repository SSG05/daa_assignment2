#include <stdio.h>
#include <stdlib.h>

void remove_first_three_lines(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening original file");
        return;
    }

    FILE *temp = fopen("temp.txt", "w");
    if (!temp) {
        perror("Error opening temporary file");
        fclose(file);
        return;
    }

    char buffer[1024];
    int line_count = 0;

    
    while (line_count < 3 && fgets(buffer, sizeof(buffer), file)) {
        line_count++;
    }

    while (fgets(buffer, sizeof(buffer), file)) {
        fputs(buffer, temp);
    }

    fclose(file);
    fclose(temp);

   
    if (remove(filename) != 0) {
        perror("Error deleting original file");
        return;
    }
    if (rename("temp.txt", filename) != 0) {
        perror("Error renaming temp file");
        return;
    }

    printf("First 3 lines removed successfully.\n");
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <filename>\n", argv[0]);
        return EXIT_FAILURE;
    }

    
    remove_first_three_lines(argv[1]);
    return EXIT_SUCCESS;
}