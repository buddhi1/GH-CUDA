// How to compile 

// g++ readShapefile.cpp -o a.out

// **************************************************

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <set>

using namespace std;

#include "../lib/point2D.h"
#include "../lib/polygon.h"
#define FILE_COUNT 14

vector<polygon> PPTmp, QQTmp;                 // two input polygons
vector<vector<double>> PPMBRs, QQMBRs;                 // two input polygons

#include <iomanip>
#include<bits/stdc++.h>
#include "../../Rtree-Topdown-version/librtree/include/index.h"

void printPolygonVector(vector<polygon>& PP) {
  int sum = 0;
  cout << PP.size() << " polygon";
  if (PP.size() > 1) 
    cout << "s";
  cout << " with \n";
  for (int i=0; i<PP.size(); i++) {
    int count = 0;
	cout << "Polygon " << i+1 << "\n";
    for (vertex* V : PP[i].vertices(ALL)){
      count++;
	  cout << setprecision (15) << V->p << endl;
	}
    cout << "my count " << count;
    if (i<PP.size()-1) 
      cout << " + \n";
    sum += count;
  }
  if (PP.size() > 1)
    cout << " = " << sum;
  cout << " vertices\n\n";
}

long double stringToDouble(string s){
	long double num=0.0;
	int precLoc=s.find(".");
	int beforePrec=s.length()-precLoc-1;
	// cout << s.length()-precLoc-1 << endl;

	// after precision
	for(int i=s.length()-1, j=0; i>precLoc; --i, ++j){
		num+=((int)s[i]-48)/pow(10, beforePrec-j);
		cout << s[i] << " - " << num << " " << ((int)s[i]-48)/ pow(10, beforePrec-j) << " ** "<< -pow(10, s.length()-1-j) << endl;
	}
	// before precision
	for(int i=precLoc-1, j=0; i>=0; --i, j++) {
		num+=((int)s[i]-48)*pow(10, j);
	}
	cout << "out " << setprecision (15) << num << endl;
	return num;
}

void loadPolygonFromShapeFile(vector<polygon>& PP, string s, int endOfFile) {
	string line;
	// cout << s << endl;
	ifstream from(s);
	// cout << "line " << line << endl;

	bool polygonReady=false;
	bool polygonStart=false;
	bool vertexFound=false;
	string polygonString="";
	string vertex="";
	string vertex2="";

	point2D v;
    polygon P;
	int count=0;

	do{
		from >> line;
		cout << line << endl;

		// check if there is comma to find vertices in the polygon
		if (polygonStart && line.find(",")!= std::string::npos) {
			vertex2=line.substr(0, line.find(","));
			vertexFound=true;
		}
		// adding end of the polygon 
		if (polygonStart && line.find("))")!= std::string::npos) {
			vertex2= line.substr(0, line.find("))"));
			polygonReady=true;
			polygonStart=false;

			PP.push_back(P);

			P = polygon(); 
			vertexFound=false;
			break;
		}
		
		if(polygonStart && !polygonReady && !vertexFound){
			vertex=line+" ";
		}
		// polygon start
		if (line.find("((")!= std::string::npos){
			vertex= line.substr(line.find("((")+2)+" ";

			polygonStart=true;
			polygonReady=false;
		} 
		
		if(vertexFound){ 
			v=point2D(atof(vertex.c_str()), atof(vertex2.c_str()));
			P.newVertex(v, true);

			vertexFound=false;
		}
		// cout << PP.size() < endOfFile << endl;
	}while((PP.size() < endOfFile || endOfFile == -1) && (!from.eof() || endOfFile != -1));

	// cout << "Num polygons " << PP.size() << endl;
	// printPolygonVector(PP);
}

// to read parks and lakes from OSM new data
void loadPolygonFromShapeFile4(vector<polygon>& PP, string s, int endOfFile, int saveOnlyId) {
	string line;
	// cout << s << endl;
	ifstream from(s);
	// cout << "line " << line << endl;

	bool polygonReady=false;
	bool polygonStart=false;
	bool vertexFound=false;
	string polygonString="";
	string vertex="";
	string vertex2="";

	point2D v;
    polygon P;
	int count=0;

	do{
		from >> line;
		// cout << line << endl;

		// check if there is comma to find vertices in the polygon
		if (polygonStart && line.find(",")!= std::string::npos) {
			vertex2=line.substr(0, line.find(","));
			vertexFound=true;
		}
		// adding end of the polygon 
		if (polygonStart && line.find("))")!= std::string::npos) {
			vertex2= line.substr(0, line.find("))"));
			polygonReady=true;
			polygonStart=false;
			if(saveOnlyId==-1 || saveOnlyId==count) PP.push_back(P);
			if(saveOnlyId==count) break;

			count++;
			P = polygon(); 
			vertexFound=false;
			// if(count++==3) break;
		}
		
		if(polygonStart && !polygonReady && !vertexFound){
			vertex=line+" ";
		}
		// polygon start
		if (line.find("((")!= std::string::npos){
			vertex= line.substr(line.find("((")+2)+" ";

			polygonStart=true;
			polygonReady=false;
		} 
		
		if(vertexFound){ 
			v=point2D(atof(vertex.c_str()), atof(vertex2.c_str()));
			P.newVertex(v, true);

			vertexFound=false;
		}
		// cout << PP.size() < endOfFile << endl;
	}while((PP.size() < endOfFile || endOfFile == -1) && (!from.eof() || endOfFile != -1));

	// cout << "Num polygons " << PP.size() << endl;
	// printPolygonVector(PP);
}

// to read parks and lakes from OSM new data and save individual polygons to files
void loadPolygonFromShapeFileNWrite4(vector<polygon>& PP, string s, int endOfFile, int saveOnlyId) {
	string line;
	// cout << s << endl;
	ifstream from(s);
	// cout << "line " << line << endl;

	bool polygonReady=false;
	bool polygonStart=false;
	bool vertexFound=false;
	string polygonString="";
	string vertex="";
	string vertex2="";

	point2D v;
    polygon P;
	int count=0;

	do{
		from >> line;
		// cout << line << endl;

		// check if there is comma to find vertices in the polygon
		if (polygonStart && line.find(",")!= std::string::npos) {
			vertex2=line.substr(0, line.find(","));
			vertexFound=true;
		}
		// adding end of the polygon 
		if (polygonStart && line.find("))")!= std::string::npos) {
			vertex2= line.substr(0, line.find("))"));
			polygonReady=true;
			polygonStart=false;
			if(saveOnlyId==-1 || saveOnlyId==count) PP.push_back(P);
			if(saveOnlyId==count) break;

			count++;
			P = polygon(); 
			vertexFound=false;
			// if(count++==3) break;
		}
		
		if(polygonStart && !polygonReady && !vertexFound){
			vertex=line+" ";
		}
		// polygon start
		if (line.find("((")!= std::string::npos){
			vertex= line.substr(line.find("((")+2)+" ";

			polygonStart=true;
			polygonReady=false;
		} 
		
		if(vertexFound){ 
			v=point2D(atof(vertex.c_str()), atof(vertex2.c_str()));
			P.newVertex(v, true);

			vertexFound=false;
		}
		// cout << PP.size() < endOfFile << endl;
	}while((PP.size() < endOfFile || endOfFile == -1) && (!from.eof() || endOfFile != -1));

	// cout << "Num polygons " << PP.size() << endl;
	// printPolygonVector(PP);
}

// for ocean dataset
void loadPolygonFromShapeFile2(vector<polygon>& PP, string s, int endOfFile) {
	string line;
	// cout << s << endl;
	ifstream from(s);
	// cout << "line " << line << endl;

	bool polygonReady=false;
	bool polygonStart=false;
	bool vertexFound=false;
	string polygonString="";
	string vertex="";
	string vertex2="", vertexNew="";

	point2D v;
    polygon P;
	int count=0;

	do{
		from >> line;
		// cout << line << endl;

		// check if there is comma to find vertices in the polygon
		if (polygonStart && line.find("),")== std::string::npos) {
			if (polygonStart && line.find(",")!= std::string::npos) {			
				vertex2=line.substr(0, line.find(","));
				vertex=vertexNew;
				vertexNew=line.substr(line.find(",")+1)+" ";
				vertexFound=true;
				// cout << line << " " << vertex2 << endl;
				// break;
			}
		}
		// adding end of the polygon 
		if (polygonStart && line.find(")))")!= std::string::npos) {
			vertex2= line.substr(0, line.find(")))"));
			polygonReady=true;
			polygonStart=false;

			PP.push_back(P);

			P = polygon(); 
			vertexFound=false;
		} else if (polygonStart && line.find(")")!= std::string::npos) {
			vertex2= line.substr(0, line.find(")"));
			polygonReady=true;
			polygonStart=false;

			PP.push_back(P);

			P = polygon(); 
			vertexFound=false;
		}
		
		// if(polygonStart && !polygonReady && !vertexFound){
		// 	vertex=line+" ";
		// }
		// polygon start
		if (line.find("(((")!= std::string::npos){
			vertexNew= line.substr(line.find("(((")+3)+" ";
			polygonStart=true;
			polygonReady=false;
			// cout << line << " " << vertex << endl;
			// break;
		} else if (line.find("(")!= std::string::npos){
			vertexNew= line.substr(line.find("(")+1)+" ";
			polygonStart=true;
			polygonReady=false;
			// cout << line << " " << vertex << endl;
			// break;
		}
		
		if(vertexFound){ 
			v=point2D(atof(vertex.c_str()), atof(vertex2.c_str()));
			P.newVertex(v, true);

			vertexFound=false;
		}
		// cout << PP.size() < endOfFile << endl;
	}while((PP.size() < endOfFile || endOfFile == -1) && (!from.eof() || endOfFile != -1));

	// cout << "Num polygons " << PP.size() << endl;
	// printPolygonVector(PP);
}

// to read ocean, land and continents from OSM new data and save individual polygons to files
void loadPolygonFromShapeFileNWrite2(vector<polygon>& PP, string s, int endOfFile, int saveOnlyId) {
	string line;
	ifstream from(s);

	bool polygonReady=false;
	bool polygonStart=false;
	bool vertexFound=false;
	string polygonString="";
	string vertex="";
	string vertex2="", vertexNew="";

	point2D v;
    polygon P;
	int count=0;

	do{
		from >> line;

		// check if there is comma to find vertices in the polygon
		if (polygonStart && line.find("),")== std::string::npos) {
			if (polygonStart && line.find(",")!= std::string::npos) {			
				vertex2=line.substr(0, line.find(","));
				vertex=vertexNew;
				vertexNew=line.substr(line.find(",")+1)+" ";
				vertexFound=true;
			}
		}

		// adding end of the polygon 
		if (polygonStart && line.find(")))")!= std::string::npos) {
			vertex2= line.substr(0, line.find(")))"));
			polygonReady=true;
			polygonStart=false;
			if(saveOnlyId==-1 || saveOnlyId==count) PP.push_back(P);
			if(saveOnlyId==count) break;

			count++;
			P = polygon(); 
			vertexFound=false;
		} else if (polygonStart && line.find(")")!= std::string::npos) {
			vertex2= line.substr(0, line.find(")"));
			polygonReady=true;
			polygonStart=false;
			if(saveOnlyId==-1 || saveOnlyId==count) PP.push_back(P);
			if(saveOnlyId==count) break;

			count++;
			P = polygon(); 
			vertexFound=false;
		}
		// polygon start
		if (line.find("(((")!= std::string::npos){
			vertexNew= line.substr(line.find("(((")+3)+" ";
			polygonStart=true;
			polygonReady=false;
		} else if (line.find("(")!= std::string::npos){
			vertexNew= line.substr(line.find("(")+1)+" ";
			polygonStart=true;
			polygonReady=false;
		}
		
		if(vertexFound){ 
			v=point2D(atof(vertex.c_str()), atof(vertex2.c_str()));
			P.newVertex(v, true);

			vertexFound=false;
		}
	}while((PP.size() < endOfFile || endOfFile == -1) && (!from.eof() || endOfFile != -1));
}

// For Continents
void loadPolygonFromShapeFile3(vector<polygon>& PP, string s, int endOfFile) {
	string line;
	// cout << s << endl;
	ifstream from(s);
	// cout << "line " << line << endl;

	bool polygonReady=false;
	bool polygonStart=false;
	bool vertexFound=false;
	string polygonString="";
	string vertex="", vertexNew="";
	string vertex2="";

	point2D v;
    polygon P;
	int count=0;

	do{
		from >> line;
		// cout << line << endl;

		// check if there is comma to find vertices in the polygon
		if (polygonStart && line.find(")),")== std::string::npos) {
			if (polygonStart && line.find(",")!= std::string::npos) {			
				vertex2=line.substr(0, line.find(","));
				vertex=vertexNew;
				vertexNew=line.substr(line.find(",")+1)+" ";
				vertexFound=true;
				// cout <<"new " << vertexNew << endl;
				// break;
			}
		}
		// adding end of the polygon 
		if (polygonStart && line.find(")))")!= std::string::npos) {
			vertex2= line.substr(0, line.find(")))"));
			polygonReady=true;
			polygonStart=false;
			// cout << "Size polygon " << P.size << endl;
			PP.push_back(P);

			P = polygon(); 
			vertexFound=false;
		} else if (polygonStart && line.find("))")!= std::string::npos) {
			vertex2= line.substr(0, line.find("))"));
			polygonReady=true;
			polygonStart=false;

			// cout << "Size* polygon " << P.size << endl;
			PP.push_back(P);
			// cout<<"*****************************************************\n\n";

			P = polygon(); 
			vertexFound=false;
			// if(PP.size()>3) break;
		}
		
		// if(polygonStart && !polygonReady && !vertexFound){
		// 	vertex=line+" ";
		// }
		// polygon start
		if (line.find("(((")!= std::string::npos){
			vertexNew=line.substr(line.find("(((")+3)+" ";
			polygonStart=true;
			polygonReady=false;
			// cout << line << " " << vertex << endl;
			// break;
		} else if (line.find("((")!= std::string::npos){
			vertexNew=line.substr(line.find("((")+2)+" ";
			polygonStart=true;
			polygonReady=false;
			// cout << line << " " << vertex << endl;
			// break;
		}
		
		if(vertexFound){ 
			v=point2D(atof(vertex.c_str()), atof(vertex2.c_str()));
			// cout<<"test v "<<v.x<<" "<<v.y<<"\n"<<endl;
			P.newVertex(v, true);

			vertexFound=false;
		}
		// cout << PP.size() < endOfFile << endl;
	}while((PP.size() < endOfFile || endOfFile == -1) && (!from.eof() || endOfFile != -1));

	// cout << "Num polygons " << PP.size() << endl;
	// printPolygonVector(PP);
}

vector<double> getMBR(polygon pol){
	vector<double> mbr;
	double maxX, minX, maxY, minY;
	maxX=pol.root->p.x;
	maxY=pol.root->p.y;
	minX=pol.root->p.x;
	minY=pol.root->p.y;

	for (vertex* V : pol.vertices(ALL)){
		if(maxX < V->p.x){
			maxX=V->p.x;
		}
		if(minX > V->p.x){
			minX=V->p.x;
		}
		if(maxY < V->p.y){
			maxY=V->p.y;
		}
		if(minY > V->p.y){
			minY=V->p.y;
		}
	//   cout << setprecision (15) << V->p << endl;
	}
	mbr.push_back(minX);
	mbr.push_back(minY);
	mbr.push_back(maxX);
	mbr.push_back(maxY);
	// uncomment to print in file format
	// printf("%f %f %f %f\n", minX, minY, maxX, maxY);
	return mbr;
}

void loadPolygonMBRs(vector<vector<double>>& PPMBRs, string s, int endOfFile) {
	string line;
	// cout << s << endl;
	ifstream from(s);
	// cout << "line " << line << endl;

	bool polygonReady=false;
	bool polygonStart=false;
	bool vertexFound=false;
	string polygonString="";
	string vertex="";
	string vertex2="";
	vector<double> mbr;

	point2D v;
    polygon P;
	int count=0;

	do{
		from >> line;

		// check if there is comma to find vertices in the polygon
		if (polygonStart && line.find(",")!= std::string::npos) {
			vertex2=line.substr(0, line.find(","));
			vertexFound=true;
		}
		// adding end of the polygon 
		if (polygonStart && line.find("))")!= std::string::npos) {
			vertex2= line.substr(0, line.find("))"));
			polygonReady=true;
			polygonStart=false;
			PPMBRs.push_back(getMBR(P));
			
			P = polygon(); 
			vertexFound=false;
		}
		
		if(polygonStart && !polygonReady && !vertexFound){
			vertex=line+" ";
		}
		// polygon start
		if (line.find("((")!= std::string::npos){
			vertex= line.substr(line.find("((")+2)+" ";

			polygonStart=true;
			polygonReady=false;
		} 
		
		if(vertexFound){ 
			v=point2D(atof(vertex.c_str()), atof(vertex2.c_str()));
			P.newVertex(v, true);

			vertexFound=false;
		}
	}while((PPMBRs.size() < endOfFile || endOfFile == -1) && (!from.eof() || endOfFile != -1));

	// cout << "Num polygons " << PPMBRs.size() << endl;
}

void loadPolygonMBRsFromVector(vector<vector<double>>& PPMBRs, vector<polygon>& PP, int endOfFile) {
	for(int i=0; i<endOfFile; ++i){
		PPMBRs.push_back(getMBR(PP[i]));
		printf("%.17g %.17g  %.17g  %.17g \n", PPMBRs[i][0], PPMBRs[i][1], PPMBRs[i][2], PPMBRs[i][3]);
		// break;
	}
}

int MySearchCallback(int id, void* arg){
	// Note: -1 to make up for the +1 when data was inserted
	printf("Hit data rect %d\n", id-1);
	return 1; // keep going
}


int main(){
	
	// open a file in write mode.
	// string fileName="parks";
	// int polyId_list[FILE_COUNT]={4, 2,1, 0};
	// int polyId_list[FILE_COUNT]={148346, 2418, 41116, 360489, 593561, 58746, 111678, 78019, 2422, 695318, 803987};
	// int polyId_list[FILE_COUNT]={684356, 10613514, 1106503, 876440, 1065174, 1243058, 943452, 844165, 169840, 321571, 140315, 173408, 1789086, 1312358, 1502030, 729699, 34579, 679995, 34622, 1512242};
	// int polyId_list[FILE_COUNT]={943452, 844165, 169840, 321571, 140315, 173408, 1789086, 1312358, 1502030, 729699, 34579, 679995, 34622, 1512242};
	int polyId_list[FILE_COUNT]={521, 1661, 1048, 1081, 1193};
	string pid;
	vector<polygon> local_poly;

	for(int fid=0; fid<FILE_COUNT; ++fid){
		
		stringstream stream;
		stream<<polyId_list[fid];
		stream>>pid;

   		ofstream outfile;
		outfile.open("continents"+pid+".txt");
		// outfile.open(fileName+"_"+pid+".txt");

		// cout<<"Writing Start: "<<fileName<<pid<<" ..."<<endl;
		cout<<"Writing Start: ..."<<endl;
		// Linux path
		// loadPolygonFromShapeFileNWrite4(local_poly, string("../../datasets/lakes_OSM_new_tiger/")+fileName, -1, polyId_list[fid]);
		// loadPolygonFromShapeFileNWrite4(local_poly, string("../../datasets/parks_OSM_new_tiger/")+fileName, -1, polyId_list[fid]);
		// loadPolygonFromShapeFile4(PPTmp, string("../../datasets/parks_OSM_new_tiger/parks"), -1);
		loadPolygonFromShapeFileNWrite2(local_poly, string("../../datasets/continents.csv"), -1, polyId_list[fid]);


		// write to file
		outfile<<local_poly[0].size<<endl;
		for (vertex* V : local_poly[0].vertices(ALL)){
			outfile<<setprecision(15)<<V->p.x<<","<<V->p.y<<endl;		
		}
		// cout<<"Writing End: "<<fileName<<polyId_list[fid]<<endl;
		cout<<"Writing End: "<<endl;
		local_poly.clear();

		outfile.close();
	}
	// loadPolygonFromShapeFile2(PPTmp, string("../../datasets/ne_10m_ocean.csv"), -1);
	// loadPolygonFromShapeFile2(QQTmp, string("../../datasets/ne_10m_land.csv"), -1);
	// loadPolygonFromShapeFile3(QQTmp, string("../../datasets/continents.csv"), -1);
	// loadPolygonFromShapeFile2(PPTmp, string("../../datasets/lake.tsv"), -1);
	// cout <<PPTmp.size()<<" polygons found"<<endl;
	// cout <<QQTmp.size()<<" polygons found"<<endl;
	// WSL / Ubuntu path
	// loadPolygonFromShapeFile(PPTmp, string("/mnt/d//Datasets/sports/sports"), -1);
	// loadPolygonFromShapeFile(PP, string("/mnt/d//Datasets/lakes/lakes"), 1);

	// print PPTmp
	// for(int ii=0; ii<PPTmp.size(); ++ii){
	// 	for (vertex* V : PPTmp[ii].vertices(ALL)){
	// 		cout<<V->p.x<<", "<<V->p.y<<" ";
	// 	}
	// 	cout<<"\n-------"<<endl;
	// }

	// load MBRs from the given polygon data
	// loadPolygonMBRs(PPMBRs, string("../../datasets/ne_10m_ocean.csv"), -1);
	// loadPolygonMBRs(PPMBRs, string("../../datasets/continents.csv"), -1);

	// load MBRs from vectors
	// loadPolygonMBRsFromVector(PPMBRs, PPTmp, PPTmp.size());
	// loadPolygonMBRsFromVector(QQMBRs, QQTmp, QQTmp.size());

	// cout<<"P Polygons" << endl;
	// // order polygon files by the number of vertices
	// vector<pair <int, int>> vect;
	// for (int i=0; i<PPTmp.size(); i++){
	// 	vect.push_back(make_pair(PPTmp[i].size, i));
	// }
	// sort(vect.rbegin(), vect.rend());
	// cout<<"============"<<endl;
	// // Printing the vector
    // for (int i=0; i<vect.size(); i++){
    // // for (int i=0; i<20; i++){
    //     cout << vect[i].first << " " << vect[i].second << endl;
    // }

	// cout<<"Q Polygons" << endl;
	// // order polygon files by the number of vertices
	// vector<pair <int, int>> vectQ;
	// for (int i=0; i<QQTmp.size(); i++){
	// 	vectQ.push_back(make_pair(QQTmp[i].size, i));
	// }
	// sort(vectQ.rbegin(), vectQ.rend());
	// // Printing the vector
    // // for (int i=0; i<vectQ.size(); i++){
    // for (int i=0; i<20; i++){
        // cout << vectQ[i].first << " " << vectQ[i].second << endl;
    // }

	// find polygon with most  number of vertices
	// int maxPolygonLoc = 0;
	// for(int i=1; i<PPTmp.size(); ++i){
	// 	if(PPTmp[i].size > PPTmp[maxPolygonLoc].size){
	// 		maxPolygonLoc=i;
	// 	}
	// }
	// cout << "Largest polygon has " << PPTmp[maxPolygonLoc].size << " Vertices at " << maxPolygonLoc << endl;

	// **use setprecision (15) to show more digits in the given double value in the output on screen
	// cout << ">> " << setprecision (15) << stod("2148.232142373") << endl;
	return 0;
}