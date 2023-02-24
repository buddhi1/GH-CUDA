vector<polygon> PPTmp, QQTmp;                 // two input polygons
vector<vector<double>> PPMBRs, QQMBRs;                 // two input polygons

#include <iomanip>
#include<bits/stdc++.h>

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

	ifstream from(s);
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
	}while((PP.size() < endOfFile || endOfFile == -1) && (!from.eof() || endOfFile != -1));
}

// to read parks and lakes from OSM new data
void loadPolygonFromShapeFile4(vector<polygon>& PP, string s, int endOfFile, int saveOnlyId) {
	string line;

	ifstream from(s);
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
	}while((PP.size() < endOfFile || endOfFile == -1) && (!from.eof() || endOfFile != -1));
}

// to read parks and lakes from OSM new data and save individual polygons to files
void loadPolygonFromShapeFileNWrite4(vector<polygon>& PP, string s, int endOfFile, int saveOnlyId) {
	string line;

	ifstream from(s);
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
	}while((PP.size() < endOfFile || endOfFile == -1) && (!from.eof() || endOfFile != -1));
}

// for ocean dataset
void loadPolygonFromShapeFile2(vector<polygon>& PP, string s, int endOfFile) {
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

	ifstream from(s);
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
		if (polygonStart && line.find(")),")== std::string::npos) {
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
			PP.push_back(P);
			P = polygon(); 
			vertexFound=false;
		} else if (polygonStart && line.find("))")!= std::string::npos) {
			vertex2= line.substr(0, line.find("))"));
			polygonReady=true;
			polygonStart=false;
			PP.push_back(P);
			P = polygon(); 
			vertexFound=false;
		}
		if (line.find("(((")!= std::string::npos){
			vertexNew=line.substr(line.find("(((")+3)+" ";
			polygonStart=true;
			polygonReady=false;
		} else if (line.find("((")!= std::string::npos){
			vertexNew=line.substr(line.find("((")+2)+" ";
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

	ifstream from(s);
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