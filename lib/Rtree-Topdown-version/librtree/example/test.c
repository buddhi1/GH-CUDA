// Compile me from ../
// gcc example/test.c libtree.a -I include/ -lm
// 
#include "../include/index.h"
// #include "index.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ENTRIES 6822
#define SEARCH_ENTRIES 1957

struct Rect rects[] = {
	{0, 0, 2, 2}, // xmin, ymin, xmax, ymax (for 2 dimensional RTree)
	{5, 5, 7, 7},
	{8, 5, 9, 6},
	{7, 1, 9, 2},
	{8, 5, 9, 6}
};
// int nrects = sizeof(rects) / sizeof(rects[0]);
struct Rect search_rect = {
	{6, 4, 10, 6}, // search will find above rects that this one overlaps
};

struct Rect rectsFromFile[ENTRIES];
struct Rect search_rectFromFile[SEARCH_ENTRIES];
int currentSearch;

int MySearchCallback(int id, void* arg) 
{
	// Note: -1 to make up for the +1 when data was inserted
	printf("%d %d\n",currentSearch, id-1);
	return 1; // keep going
}


// read csv file into array
void readShapeFileToArray(struct Rect arr[], char *filename){
	int i=0, j=0, end=1;
	char *line=NULL;
	char *ptr;

	FILE* fp = fopen(filename, "r");

	if (!fp)
		printf("Can't open file\n");

	else {
		char buffer[1024];

		// while (fgets(buffer, 64, fp)) {
		while (end == 1) {				
			// read into array
			for(j=0; j<4 && !feof(fp); ++j){
   				// printf("** %s\n", buffer);
				fscanf(fp, "%s", buffer);
				if(feof(fp)){
					end=-1;
					break;
				}
   				// printf("** %f %s [%d] ", strtod(buffer, &ptr), buffer, j);
				arr[i].boundary[j]=strtod(buffer, &ptr);
			}
			// printf("\n");
			if(feof(fp) && end != -1) {
				printf("ERROR!!! Endf of file reached before reading 4 values of MBR. Missing values in line %d\n", i);
				break;
			}

			// printf("file buffer %s", buffer);
			// line=strstr(buffer, "((");
			// if(line != NULL) {
   			// 	printf("** %s: %s\n", line+2, buffer);
				
			// }

			// line=strstr(buffer, "))");
			// if(line != NULL) {
   			// 	printf("** %s: %s %d\n", line, buffer, strlen(buffer));
			// 	buffer[strlen(buffer)-2] = '\0';
   			// 	printf("** %s: %s %d\n", line, buffer, strlen(buffer));
				
			// }
			i++;
			// if(i>=9)
			// 	break;
		}
		printf("array %d\n", i-1);
		fclose(fp);
	}
	
}


void main()
{
	struct Node* root = RTreeNewIndex();
	int i, nhits;

	readShapeFileToArray(rectsFromFile, "../../2019-Foster/polyclip2/ne_10m_ocean_MBRs.txt");
	readShapeFileToArray(search_rectFromFile, "../../2019-Foster/polyclip2/continents_MBRs.txt");
	// int nrects = sizeof(rectsFromFile) / sizeof(rectsFromFile[0].boundary[0]);
	printf("nrects = %d\n", ENTRIES);
	/*
	 * Insert all the data rects.
	 * Notes about the arguments:
	 * parameter 1 is the rect being inserted,
	 * parameter 2 is its ID. NOTE: *** ID MUST NEVER BE ZERO ***, hence the +1,
	 * parameter 3 is the root of the tree. Note: its address is passed
	 * because it can change as a result of this call, therefore no other parts
	 * of this code should stash its address since it could change undernieth.
	 * parameter 4 is always zero which means to add from the root.
	 */
	// for(i=0; i<ENTRIES; ++i){
	// 	printf("%f ", search_rectFromFile[i].boundary[0]);
	// 	printf("%f ", search_rectFromFile[i].boundary[1]);
	// 	printf("%f ", search_rectFromFile[i].boundary[2]);
	// 	printf("%f ", search_rectFromFile[i].boundary[3]);
	// 	printf("\n");
	// }

	for(i=0; i<ENTRIES; i++)
		RTreeInsertRect(&rectsFromFile[i], i+1, &root, 0); // i+1 is rect ID. Note: root can change
	
	for(i=0; i<SEARCH_ENTRIES; i++){
		nhits=0;
		currentSearch=i;
		nhits = RTreeSearch(root, &search_rectFromFile[i], MySearchCallback, 0);
		// printf("Search resulted in %d hits\n", nhits);
	}
}