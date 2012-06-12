#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
FILE *fp;
char buffer [1000];
int count, i;
fpos_t pos;
if(argc <= 2) return 0;
if ((fp = fopen(argv[1], "r")) == NULL) {
	fprintf(stderr, "Can't open %s\n", argv[1]);
	exit(1);
} 
else {
	fseek(fp, 0, SEEK_END);
	count = atoi (argv[2]);
	
	for (i=0;i<count;i++) {
		ungetc(i, fp);
	}
	fgetpos(fp, &pos);
	fsetpos(fp, &pos);
	fgets (buffer, count, fp);
	printf("%s Characters Read from file %s : \n%s\n", argv[2], argv[1], buffer);

	fclose(fp);
}
}
