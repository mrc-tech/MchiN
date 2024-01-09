#include <iostream>
#include <fstream>
#include <string>

using namespace std;



int main(int argc,char** argv)
{
	ifstream input;
	ofstream output("output.txt");
	
	input.open(argv[1]);
	
	string str;
	while(getline(input,str)){
		string temp;
		for(int i=0;i<str.length();i++){
			if(str[i]=='à') temp = temp + "a\'";
			else if(str[i]=='è') temp = temp + "e\'";
			else if(str[i]=='é') temp = temp + "e\'";
			else if(str[i]=='ì') temp = temp + "i\'";
			else if(str[i]=='ò') temp = temp + "o\'";
			else if(str[i]=='ù') temp = temp + "u\'";
			else temp.push_back(str[i]);
		}
		cout << temp << endl;
		output << temp << endl;
	}
	
	input.close();
	output.close();
	
	
	return 0;
}
