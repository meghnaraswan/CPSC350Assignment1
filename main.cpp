//  Meghna Raswan
//  2337415
//  raswan@chapman.edu
//  CPSC 350
//  Assignment 1
//  main.cpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include<cmath>
#include <math.h>
#include <algorithm>

using namespace std;

class MyDNA {
public:
    //member variables

    //member variables for the character count of the file (or the sum of each individual nucleotides) anf the number of lines
    int m_charCount = 0;
    int m_lineCount = 0;
    
    //member variables for each nucleotide
    int m_C = 0;
    int m_A = 0;
    int m_T = 0;
    int m_G = 0;
    
    //member variables for the mean, variance, and standard deviation
    double m_mean = 0.0;
    double m_variance = 0.0;
    double m_stddev = 0.0;
    
    //member variables for the realtive probability of each nucleotide
    double m_CProb;
    double m_AProb;
    double m_TProb;
    double m_GProb;
    
    //member variables for each bigram
    int m_AA = 0;
    int m_AC = 0;
    int m_AT = 0;
    int m_AG = 0;
    int m_CA = 0;
    int m_CC = 0;
    int m_CT = 0;
    int m_CG = 0;
    int m_TA = 0;
    int m_TC = 0;
    int m_TT = 0;
    int m_TG = 0;
    int m_GA = 0;
    int m_GC = 0;
    int m_GT = 0;
    int m_GG = 0;
    
    //member variable for the total amount of bigrams
    int m_TotalBigrams = 0;
    
    //member variables for the realtive probability of each nucleotide bigram
    double m_AAProb;
    double m_ACProb;
    double m_ATProb;
    double m_AGProb;
    double m_CAProb;
    double m_CCProb;
    double m_CTProb;
    double m_CGProb;
    double m_TAProb;
    double m_TCProb;
    double m_TTProb;
    double m_TGProb;
    double m_GAProb;
    double m_GCProb;
    double m_GTProb;
    double m_GGProb;
    
    //member variables for input and output files
    string m_filename;
    string m_outputfile;
    
public:
    //public constructor for input and output files
    MyDNA(string filename, string outputfile) {
        m_filename = filename;
        m_outputfile = outputfile;
    }
    
    //function for calculating the mean, character counts, and number of each nucleotide
    void calcMeanAndCounts(){
        ifstream ifs (m_filename); //constructs ifstream object and opens user input file
        if (ifs.is_open()) { //shows if the stream is associated with the file
            char c = ifs.get(); //extracts character from the stream
            while (ifs.good()) { //checks the state of the stream (so checks to see if the stream is good)
                switch(c){ //switch statement for counting the number of nucleotides and each individual nucleotide depending on the character value in the file
                    case 'A':
                    case 'a':
                        //if character value is 'A' or 'a', then this increments the member variable and character count by 1
                        ++m_A;
                        ++m_charCount;
                        break;
                    case 'C':
                    case 'c':
                        //if character value is 'C' or 'c', then this increments the member variable and character count by 1
                        ++m_C;
                        ++m_charCount;
                        break;
                    case 'T':
                    case 't':
                        //if character value is 'T' or 't', then this increments the member variable and character count by 1
                        ++m_T;
                        ++m_charCount;
                        break;
                    case 'G':
                    case 'g':
                        //if character value is 'G' or 'g', then this increments the member variable and character count by 1
                        ++m_G;
                        ++m_charCount;
                        break;
                    case '\n':
                        //if character value is '\n' (which means new line), then this increments the line count by 1
                        ++m_lineCount;
                    default:
                        break;
                }
                c = ifs.get(); //extracts character from the stream
            }
            ifs.close(); //closes file
        } else { //if file can't open, there will be an error message
            std::cout << "Error opening file: " << m_filename << endl;
        }
        m_mean = double(m_charCount)/double(m_lineCount); //calculate mean by dividing character count by line count
    }

    //function for the mean, whcih returns the mean member variable
    double getMean() {
        return m_mean;
    }
    
    //function for calculutating variance and standard deviation
    void calcVarianceStdDev(){
        ifstream ifs (m_filename);
        int lineLength = 0; //initialize the length of line to 0
        double sum_of_sqr = 0.0; //initialize the length of sum of squares to 0.0
        if (ifs.is_open()) {
            char c = ifs.get();
            while (ifs.good()) {
                switch(c){ //switch statement for adding counting the length of the line by adding each nucleotide in each line to the lineLength variable
                    case 'A':
                    case 'a':
                    case 'C':
                    case 'c':
                    case 'T':
                    case 't':
                    case 'G':
                    case 'g':
                        ++lineLength;
                        break;
                    case '\n': //in case of a new line, we will calculate each expression in the variance, ultimately resulting in the calculation of the variance after the last line (we will only be calculating for the numerator here)
                        sum_of_sqr += pow(lineLength - getMean(),2); //numerator is the difference between the length of each line and the mean squared
                        lineLength = 0; //reassign the length of the line to 0 after computing for each expression, which will continutally add to the sum of the squares equation until the sequence is over
                    default:
                        break;
                }
                c = ifs.get();
            }
            ifs.close();
        } else {
            std::cout << "Error opening file";
        }
        m_variance = sum_of_sqr / (m_lineCount); //variance is the sum of the squares divided by the number of lines
        m_stddev = sqrt(m_variance); //standard deviation is the square root of the variance
        ProbabilityOfEachNucleotide(); //calling ProbabilityOfEachNucleotide() function
    }

    //function for calculating the relative probability of each nucleotide
    void ProbabilityOfEachNucleotide() { //to calculate the relative proability, each nucleotide occurence will divide by the sum of all of the nucleotides
        m_CProb = (double)m_C / (double)m_charCount;
        m_AProb = (double)m_A / (double)m_charCount;
        m_TProb = (double)m_T / (double)m_charCount;
        m_GProb = (double)m_G / (double)m_charCount;
    }
    
    //function for the bigrams in the file
    void Bigrams(string sBigram){
        //bigram sting will be the parameter for this function
        //if the string is a bigram that matches up with any of these bigrams listed below, the member variables for each individual bigram and the total bigrams will increase by 1
        if (sBigram == "AA"){
            ++m_AA;
            ++m_TotalBigrams;
        }
        else if (sBigram == "AC"){
            ++m_AC;
            ++m_TotalBigrams;
        }
        else if (sBigram == "AT"){
            ++m_AT;
            ++m_TotalBigrams;
        }
        else if (sBigram == "AG"){
            ++m_AG;
            ++m_TotalBigrams;
        }
        else if (sBigram == "CA"){
            ++m_CA;
            ++m_TotalBigrams;
        }
        else if (sBigram == "CC"){
            ++m_CC;
            ++m_TotalBigrams;
        }
        else if (sBigram == "CT"){
            ++m_CT;
            ++m_TotalBigrams;
        }
        else if (sBigram == "CG"){
            ++m_CG;
            ++m_TotalBigrams;
        }
        else if (sBigram == "TA"){
            ++m_TA;
            ++m_TotalBigrams;
        }
        else if (sBigram == "TC"){
            ++m_TC;
            ++m_TotalBigrams;
        }
        else if (sBigram == "TT"){
            ++m_TT;
            ++m_TotalBigrams;
        }
        else if (sBigram == "TG"){
            ++m_TG;
            ++m_TotalBigrams;
        }
        else if (sBigram == "GA"){
            ++m_GA;
            ++m_TotalBigrams;
        }
        else if (sBigram == "GC"){
            ++m_GC;
            ++m_TotalBigrams;
        }
        else if (sBigram == "GT"){
            ++m_GT;
            ++m_TotalBigrams;
        }
        else if (sBigram == "GG"){
            ++m_GG;
            ++m_TotalBigrams;
        }
    }

    //function for counting the bigrams
    void CountBigram() {
        string s; //creating an s string variable
        ifstream ifs (m_filename);
        if (ifs.is_open()) {
            getline(ifs, s); //reads each line in the file
            s.erase(remove(s.begin(), s.end(), '\r'), s.end()); //this removes '\r' from the file; .erase() removes '\r' reanging from the beginning to the end
            while (ifs.good()) {
                if (s.length() % 2 == 0){ //if the length of the string is even
                    for (int i = 0; i < s.length(); ++i, ++i){ //iterate twice since we will be looking at 2 strings at a time
                        string sBigram = ""; //assign sBigram as an empty string
                        sBigram.push_back(s[i]); //push_back will add new element at the end of the string, or verctor; first index of string will be even
                        sBigram.push_back(s[i + 1]); //second index of string will be odd
                        Bigrams(sBigram); //calling the Brigrams() function, passing sBigram as a parameter
                    }
                }
                else { //in case the length of the string is odd; same ules still apply
                    for (int i = 0; i < s.length(); ++i){
                        string sBigram = "";
                        sBigram.push_back(s[i]);
                        sBigram.push_back(s[i + 1]);
                        Bigrams(sBigram);
                    }
                }
                getline(ifs, s);
                s.erase( remove(s.begin(), s.end(), '\r'), s.end() ); //removes '\r' from the file
            }
            ifs.close();
        } else {
            std::cout << "Error opening file";
        }
    }
    
    //function for calculating the relative probability of each nucleotide bigram
    void BigramRelativeProability() {
        //to calculate the relative proability, each bigram occurence will divide by the sum of all of the bigrams
        m_AAProb = (double)m_AA / (double)m_TotalBigrams;
        m_ACProb = (double)m_AC / (double)m_TotalBigrams;
        m_ATProb = (double)m_AT / (double)m_TotalBigrams;
        m_AGProb = (double)m_AG / (double)m_TotalBigrams;
        m_CAProb = (double)m_CA / (double)m_TotalBigrams;
        m_CCProb = (double)m_CC / (double)m_TotalBigrams;
        m_CTProb = (double)m_CT / (double)m_TotalBigrams;
        m_CGProb = (double)m_CG / (double)m_TotalBigrams;
        m_TAProb = (double)m_TA / (double)m_TotalBigrams;
        m_TCProb = (double)m_TC / (double)m_TotalBigrams;
        m_TTProb = (double)m_TT / (double)m_TotalBigrams;
        m_TGProb = (double)m_TG / (double)m_TotalBigrams;
        m_GAProb = (double)m_GA / (double)m_TotalBigrams;
        m_GCProb = (double)m_GC / (double)m_TotalBigrams;
        m_GTProb = (double)m_GT / (double)m_TotalBigrams;
        m_GGProb = (double)m_GG / (double)m_TotalBigrams;
    }
    
    //function for converting all of the data into a string
    string ToString(){
        string ret_Str = "Line Count: " + to_string(m_lineCount) + "\n";
        ret_Str +=  "Sum: " + to_string(m_charCount) + "\n";
        ret_Str +=  "C: " + to_string(m_C) + "\n";
        ret_Str +=  "G: " + to_string(m_G) + "\n";
        ret_Str +=  "A: " + to_string(m_A) + "\n";
        ret_Str +=  "T: " + to_string(m_T) + "\n";
        ret_Str +=  "Mean: " + to_string(m_mean) + "\n";
        ret_Str +=  "Variance: " + to_string(m_variance) + "\n";
        ret_Str +=  "Standard Deviation: " + to_string(m_stddev) + "\n";
        ret_Str +=  "Relative probability for each nucleotide: \n";
        ret_Str +=  "    A: " + to_string(m_AProb) + "\n";
        ret_Str +=  "    C: " + to_string(m_CProb) + "\n";
        ret_Str +=  "    T: " + to_string(m_TProb) + "\n";
        ret_Str +=  "    G: " + to_string(m_GProb) + "\n";
        ret_Str +=  "Relative probability of each nucleotide bigram: \n";
        ret_Str +=  "    AA: " + to_string(m_AAProb) + "\n";
        ret_Str +=  "    AC: " + to_string(m_ACProb) + "\n";
        ret_Str +=  "    AT: " + to_string(m_ATProb) + "\n";
        ret_Str +=  "    AG: " + to_string(m_AGProb) + "\n";
        ret_Str +=  "    CA: " + to_string(m_CAProb) + "\n";
        ret_Str +=  "    CC: " + to_string(m_CCProb) + "\n";
        ret_Str +=  "    CT: " + to_string(m_CTProb) + "\n";
        ret_Str +=  "    CG: " + to_string(m_CGProb) + "\n";
        ret_Str +=  "    TA: " + to_string(m_TAProb) + "\n";
        ret_Str +=  "    TC: " + to_string(m_TCProb) + "\n";
        ret_Str +=  "    TT: " + to_string(m_TTProb) + "\n";
        ret_Str +=  "    TG: " + to_string(m_TGProb) + "\n";
        ret_Str +=  "    GA: " + to_string(m_GAProb) + "\n";
        ret_Str +=  "    GC: " + to_string(m_GCProb) + "\n";
        ret_Str +=  "    GT: " + to_string(m_GTProb) + "\n";
        ret_Str +=  "    GG: " + to_string(m_GGProb) + "\n";
        return ret_Str;
    }
    
    //print statement for all of the information
    void Print(){
        cout << "Line Count: " << m_lineCount << endl;
        cout << "Sum: " << m_charCount << endl;
        cout << "C: " << m_C << endl;
        cout << "G: " << m_G << endl;
        cout << "A: " << m_A << endl;
        cout << "T: " << m_T << endl;
        cout << "Mean: " << m_mean << endl;
        cout << "Variance: " << m_variance << endl;
        cout << "Standard Deviation: " << m_stddev << endl;
        cout << "Relative probability for each nucleotide: " << endl;
        cout << "    A: " << m_AProb << endl;
        cout << "    C: " << m_CProb << endl;
        cout << "    T: " << m_TProb << endl;
        cout << "    G: " << m_GProb << endl;
        cout << "Relative probability of each nucleotide bigram:" << endl;
        cout << "    AA: " << m_AAProb << endl;
        cout << "    AC: " << m_ACProb << endl;
        cout << "    AT: " << m_ATProb << endl;
        cout << "    AG: " << m_AGProb << endl;
        cout << "    CA: " << m_CAProb << endl;
        cout << "    CC: " << m_CCProb << endl;
        cout << "    CT: " << m_CTProb << endl;
        cout << "    CG: " << m_CGProb << endl;
        cout << "    TA: " << m_TAProb << endl;
        cout << "    TC: " << m_TCProb << endl;
        cout << "    TT: " << m_TTProb << endl;
        cout << "    TG: " << m_TGProb << endl;
        cout << "    GA: " << m_GAProb << endl;
        cout << "    GC: " << m_GCProb << endl;
        cout << "    GT: " << m_GTProb << endl;
        cout << "    GG: " << m_GGProb << endl;
    }
    
    //function for the gaussian distribution; returns integer value
    int GaussianDistribution() {
        double randomVariableC; //first random variable C
        double randomVariableD; //second random variable D
        int lengthOfSequence; //length of the sequence variable
        double a; //variable a
        double b; //variable b
        a = ((double) rand() / (RAND_MAX)); //random number divides the max number which will always assign variable "a" a number between 0 and 1
        b = ((double) rand() / (RAND_MAX)); //random number divides the max number which will always assign variable "b" a number between 0 and 1
        randomVariableC = sqrt(-2 * log10(a)) * cos(2 * M_PI * b); //Box-Muller transport that will assign C a number followed by this equation
        randomVariableD = m_stddev * randomVariableC + m_mean; //random number D is the product of the tandard deviation and random number C plus the mean
        lengthOfSequence = int(randomVariableD); //assign the length of the sequence as the integer value of random number D
        return lengthOfSequence; //return length of sequence
    }
    
    //function for printing all of the information into an output file
    void PrintingInFile() {
        int lineCount = 0; //initialize line count to 0
        ofstream out(m_outputfile); //output stream class to operate on files
        
        //outputs all of this info into a file
        out << "Meghna Raswan\n";
        out << "2337415\n";
        out << "raswan@chapman.edu\n";
        out << "CPSC 350\n";
        out << "Assignment 1\n";
        
        do {
            int myLineLength = GaussianDistribution(); //call GaussianDistribution() to find length of the line of the sequence (using the return value of GaussianDistribution() function to find the length of the line)
            out << myLineLength << " - " ;
            
            //to find the length of each nucleotide, find product of the length of the line and the relative proability fo each nucleotide
            //use forloops for each nucleotide to find ech of their lengths and add them to the file on each line
            int myALength = round(myLineLength * m_AProb);
            for (int i = 0; i < myALength; ++i){
                out.put('A');
            }
            int myTLength = round(myLineLength * m_TProb);
            for (int i = 0; i < myTLength; ++i){
                out.put('T');
            }
            int myCLength = round(myLineLength * m_CProb);
            for (int i = 0; i < myCLength; ++i){
                out.put('C');
            }
            int myGLength = round(myLineLength * m_GProb);
            for (int i = 0; i < myGLength; ++i){
                out.put('G');
            }
            ++ lineCount; //increment th e line count by 1 after the length of the line if finished computing
            out << "\n";
        } while (lineCount < 1000); //this will pring 1000 lines
        out << this->ToString(); //call the ToString() method with a pointer
        out.close(); //close the output file
    }
};

//main function
int main(int argc, char *argv[]) {
    string userFileName; //user file name variable
    if (argc < 2) { //if the arguments on the command line are less than 2 (meaning the user did not add the file name), then return an error
        cout << "usage: ComputeDNASeq <INPUT_FILE_NAME>";
        return 0; //exit
    } else { //else, userFileName is the commandline user input
        userFileName = argv[1];
    }
    cout << "Input File :" << userFileName <<"\n";
    while (true){ //allows the program to keep running until the user exits the program
        string userOutputFile = userFileName + ".out.txt"; //output file
        MyDNA myDna(userFileName, userOutputFile); //creating opject with inout and output file parameters
        myDna.calcMeanAndCounts(); //calling calcMeanAndCounts() function
        myDna.calcVarianceStdDev(); //calling calcVarianceStdDev() function
        myDna.CountBigram(); //calling CountBigram() function
        myDna.BigramRelativeProability(); //calling BigramRelativeProability() function
        
        myDna.Print(); //calling Print() function
        myDna.PrintingInFile(); //calling PrintingInFile() function
        
        string option; //option for if user wants to continue with examining more inout files or if user wants to exit program
        cout << "Would you still like to continue to compute the DNA Sequence info? Press 'Y' to continue or press 'N' to quit:" << endl;
        cin >> option; //user inout for option
        if (option == "Y" || option == "y" ) { //if user enters "Y" or "y", then the program continues
            cout << "Input File: ";
            cin >> userFileName; //user inout
            continue;
        } else { //else the program exits
            cout << "Bye! \n";
            return 0;
        }
    }
    return 0;
}

