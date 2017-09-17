#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#define CSVFILEPATH "./csvDump/"
#define TXTFILEPATH "./txtDump/"

typedef std::vector<std::vector<int>> arrayOfInt2D;
typedef std::vector<std::vector<std::vector<int>>> arrayOfInt3D;
typedef std::vector<std::vector<std::vector<std::vector<int>>>> arrayOfInt4D;
typedef std::vector<std::vector<double>> arrayOfDouble2D;
typedef std::vector<std::vector<std::vector<double>>> arrayOfDouble3D;
typedef std::vector<std::vector<std::vector<std::vector<double>>>> arrayOfDouble4D;

void fun();

arrayOfInt2D getArrayOfInt2D(int n, int m, int val = 0);
arrayOfInt3D getarrayOfInt3D(int n, int m, int p, int val = 0);
arrayOfInt4D getArrayOfInt4D(int n, int m, int p, int q, int val = 0);
arrayOfDouble2D getArrayOfDouble2D(int n, int m, double val = 0.0);
arrayOfDouble3D getArrayOfDouble3D(int n, int m, int p, double val = 0.0);
arrayOfDouble4D getArrayOfDouble4D(int n, int m, int p, int q, double val = 0.0);

double getMinimumOf2DArray(arrayOfDouble2D array2D);
double getMinimumOf3DArray(arrayOfDouble3D array3D, int& c);

double getMaximumOfArray(std::vector<double> vec);
double getMaximumOf2DArray(arrayOfDouble2D array2D);

int getCountOfNegativeIn3DArray(arrayOfDouble3D array3D);

std::vector<double> linearize3DVector(arrayOfDouble3D array3D);

std::string moreSigs(double d, int prec); //return string of 'd' with ''prec' sig digits: trailing zeros removed

std::vector<std::string> listFiles(std::string path, std::string ext);

template <typename T>
void dumpData(T data, std::string varName)
{
    std::string path = TXTFILEPATH;
    std::string fileName = varName + ".txt";
    std::ofstream myFile;
    myFile.open((path + fileName).c_str());
    if (!myFile.is_open())
    {
        std::cout << "Unable to open file to dump data" << std::endl;
        return;
    }
    for (size_t i = 0; i < data.size(); i++)
        myFile << data[i] << std::endl;

    myFile.close();
}

template <typename T>
void dump2DData(T data, std::string varName)
{
    std::string path = TXTFILEPATH;
    std::string fileName = varName + ".txt";
    std::ofstream myFile;
    myFile.open((path + fileName).c_str());
    if (!myFile.is_open())
    {
        std::cout << "Unable to open file to dump data" << std::endl;
        return;
    }

    for (size_t i = 0; i < data.size(); i++)
    {
        myFile << "i = " << i << std::endl;
        for (size_t j = 0; j < data[i].size(); j++)
            myFile << data[i][j] << '\t';
        myFile << std::endl
               << std::endl;
    }
    myFile.close();
}

template <typename T>
void dump3DData(T data, std::string varName)
{
    std::string path = TXTFILEPATH;
    std::string fileName = varName + ".txt";
    std::ofstream myFile;
    myFile.open((path + fileName).c_str());
    if (!myFile.is_open())
    {
        std::cout << "Unable to open file to dump data" << std::endl;
        return;
    }

    size_t thirdDim = data[0][0].size();
    for (size_t d3 = 0; d3 < thirdDim; d3++)
    {
        myFile << varName.c_str() << "(:,:" << d3 << ")" << std::endl;
        for (size_t d1 = 0; d1 < data.size(); d1++)
        {
            for (size_t d2 = 0; d2 < data[d1].size(); d2++)
                myFile << data[d1][d2][d3] << '\t';
            myFile << std::endl;
        }
        myFile << std::endl
               << std::endl;
    }
    myFile.close();
}

template <typename T>
void dump4DData(T data, std::string varName)
{
    std::string path = TXTFILEPATH;
    std::string fileName = varName + ".txt";
    std::ofstream myFile;
    myFile.open((path + fileName).c_str());
    if (!myFile.is_open())
    {
        std::cout << "Unable to open file to dump data" << std::endl;
        return;
    }

    size_t thirdDim = data[0][0].size();
    size_t forthDim = data[0][0][0].size();
    for (size_t d4 = 0; d4 < forthDim; d4++)
    {
        for (size_t d3 = 0; d3 < thirdDim; d3++)
        {
            myFile << varName.c_str() << "(:,:" << d3 << "," << d4 << ")" << std::endl;
            for (size_t d1 = 0; d1 < data.size(); d1++)
            {
                for (size_t d2 = 0; d2 < data[d1].size(); d2++)
                    myFile << data[d1][d2][d3][d4] << '\t';
                myFile << std::endl;
            }
            myFile << std::endl
                   << std::endl;
        }
    }
    myFile.close();
}

//Dumping data as column

template <typename T>
void dumpCSV(T data, std::string varName)
{
    std::string path = CSVFILEPATH;
    std::string fileName = varName + ".csv";
    std::ofstream myFile;
    myFile.open((path + fileName).c_str());
    if (!myFile.is_open())
    {
        std::cout << "Unable to open file to dump data" << std::endl;
        return;
    }
    myFile << "dim_1"
           << ","
           << "value" << std::endl;
    for (size_t i = 0; i < data.size(); i++)
    {
        if (data[i] < 1.0e-16)
            myFile << i + 1 << "," << data[i] << std::endl;
        else
            myFile << i + 1 << "," << moreSigs(data[i], 16) << std::endl;
    }
    myFile.close();
}

template <typename T>
void dump2DCSV(T data, std::string varName)
{
    std::string path = CSVFILEPATH;
    std::string fileName = varName + ".csv";
    std::ofstream myFile;
    myFile.open((path + fileName).c_str());
    if (!myFile.is_open())
    {
        std::cout << "Unable to open file to dump data" << std::endl;
        return;
    }
    myFile << "dim_1"
           << ","
           << "dim_2"
           << ","
           << "value" << std::endl;
    for (size_t i = 0; i < data.size(); i++)
        for (size_t j = 0; j < data[i].size(); j++)
        {
            if (data[i][j] < 1.0e-16)
                myFile << i + 1 << "," << j + 1 << "," << data[i][j] << std::endl;
            else
                myFile << i + 1 << "," << j + 1 << "," << moreSigs(data[i][j], 16) << std::endl;
        }
    myFile.close();
}

template <typename T>
void dump3DCSV(T data, std::string varName)
{
    std::string path = CSVFILEPATH;
    std::string fileName = varName + ".csv";
    std::ofstream myFile;
    myFile.open((path + fileName).c_str());
    if (!myFile.is_open())
    {
        std::cout << "Unable to open file to dump data" << std::endl;
        return;
    }
    myFile << "dim_1"
           << ","
           << "dim_2"
           << ","
           << "dim_3"
           << ","
           << "value" << std::endl;
    for (size_t d1 = 0; d1 < data.size(); d1++)
        for (size_t d2 = 0; d2 < data[d1].size(); d2++)
            for (size_t d3 = 0; d3 < data[d1][d2].size(); d3++)
            {
                if (data[d1][d2][d3] < 1.0e-16)
                    myFile << d1 + 1 << "," << d2 + 1 << "," << d3 + 1 << "," << data[d1][d2][d3] << std::endl;
                else
                    myFile << d1 + 1 << "," << d2 + 1 << "," << d3 + 1 << "," << moreSigs(data[d1][d2][d3], 16) << std::endl;
            }
    myFile.close();
}

template <typename T>
void dump4DCSV(T data, std::string varName)
{
    std::string path = CSVFILEPATH;
    std::string fileName = varName + ".csv";
    std::ofstream myFile;
    myFile.open((path + fileName).c_str());
    if (!myFile.is_open())
    {
        std::cout << "Unable to open file to dump data" << std::endl;
        return;
    }
    myFile << "dim_1"
           << ","
           << "dim_2"
           << ","
           << "dim_3"
           << ","
           << "dim_4"
           << ","
           << "value" << std::endl;
    for (size_t d1 = 0; d1 < data.size(); d1++)
        for (size_t d2 = 0; d2 < data[d1].size(); d2++)
            for (size_t d3 = 0; d3 < data[d1][d2].size(); d3++)
                for (size_t d4 = 0; d4 < data[d1][d2][d3].size(); d4++)
                {
                    if (data[d1][d2][d3][d4] < 1.0e-16)
                        myFile << d1 + 1 << "," << d2 + 1 << "," << d3 + 1 << "," << d4 + 1 << "," << data[d1][d2][d3][d4] << std::endl;
                    else
                        myFile << d1 + 1 << "," << d2 + 1 << "," << d3 + 1 << "," << d4 + 1 << "," << moreSigs(data[d1][d2][d3][d4], 16) << std::endl;
                }
    myFile.close();
}

template <typename T>
void dumpTestCSV(T data, std::string varName, double aggKerVal = 0.0, double aggKerConst = 0.0, double brkKerVal = 0.0, double brkKerConst = 0.0)
{
    std::string path = CSVFILEPATH;
    std::string fileName = varName + ".csv";
    std::ofstream myFile;
    myFile.open((path + fileName).c_str());
    if (!myFile.is_open())
    {
        std::cout << "Unable to open file to dump data" << std::endl;
        return;
    }
    myFile << ""
           << ","
           << ""
           << ","
           << ""
           << "," << aggKerVal << std::endl;
    myFile << ""
           << ","
           << ""
           << ","
           << ""
           << "," << aggKerConst << std::endl;
    myFile << ""
           << ","
           << ""
           << ","
           << ""
           << "," << brkKerVal << std::endl;
    myFile << ""
           << ","
           << ""
           << ","
           << ""
           << "," << brkKerConst << std::endl;
    myFile << "dim_1"
           << ","
           << "dim_2"
           << ","
           << "dim_3"
           << ","
           << "value" << std::endl;
    for (size_t d1 = 0; d1 < data.size(); d1++)
        for (size_t d2 = 0; d2 < data[d1].size(); d2++)
            for (size_t d3 = 0; d3 < data[d1][d2].size(); d3++)
            {
                if (data[d1][d2][d3] < 1.0e-16)
                    myFile << d1 + 1 << "," << d2 + 1 << "," << d3 + 1 << "," << data[d1][d2][d3] << std::endl;
                else
                    myFile << d1 + 1 << "," << d2 + 1 << "," << d3 + 1 << "," << moreSigs(data[d1][d2][d3], 16) << std::endl;
            }
    myFile.close();
}

template <typename T1, typename T2> //T1 for time & T2 for d10, d50 & d90 data
void dumpDiaCSV(T1 data1, T2 data2, std::string varName)
{
    std::string path = CSVFILEPATH;
    std::string fileName = varName + ".csv";
    std::ofstream myFile;
    myFile.open((path + fileName).c_str());
    if (!myFile.is_open())
    {
        std::cout << "Unable to open file to dump data" << std::endl;
        return;
    }
    myFile << "Time Index"
           << ","
           << "Time";
    for (size_t c = 0; c < data2[0].size(); c++)
        myFile << ","
               << "C" << c + 1;
    myFile << std::endl;

    for (size_t t = 0; t < data1.size(); t++)
    {
        myFile << t + 1 << "," << data1[t];
        for (size_t c = 0; c < data2[t].size(); c++)
            myFile << "," << moreSigs(data2[t][c], 16);
        myFile << std::endl;
    }
    myFile.close();
}

template <typename T>
void dump2DCSV4Matlab(T data, std::string varName)
{
    std::string path = CSVFILEPATH;
    std::string fileName = varName + ".csv";
    std::ofstream myFile;
    myFile.open((path + fileName).c_str());
    if (!myFile.is_open())
    {
        std::cout << "Unable to open file to dump data" << std::endl;
        return;
    }

    for (size_t i = 0; i < data.size(); i++)
    {

        for (size_t j = 0; j < data[i].size(); j++)
        {
            //if (data[i][j] < 1.0e-16)
            myFile << data[i][j];
            if (j < data[i].size() - 1)
                myFile << ",";
            //else
            //myFile << moreSigs(data[i][j], 16) <<",";
        }
        myFile << std::endl;
    }
    myFile.close();
}

#endif // UTILITY_H
