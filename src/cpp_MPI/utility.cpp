#include <iostream>
#include <cmath>
#include <float.h>
#include <sstream>
#include <iomanip>
#include <dirent.h>
#include <cstring>
#include "utility.h"

using namespace std;

void fun()
{
    cout << "utility::fun is called" << endl;
}

arrayOfInt2D getArrayOfInt2D(int n, int m, int val)
{
    arrayOfInt2D arrayInt2D(n, vector<int>(m, val));
    return arrayInt2D;
}

arrayOfInt3D getarrayOfInt3D(int n, int m, int p, int val)
{
    arrayOfInt3D arrayInt3D(n, vector<vector<int>>(m, vector<int>(p, val)));
    return arrayInt3D;
}

arrayOfInt4D getArrayOfInt4D(int n, int m, int p, int q, int val)
{
    arrayOfInt4D arrayInt4D(n, vector<vector<vector<int>>>(m, vector<vector<int>>(p, vector<int>(q, val))));
    return arrayInt4D;
}

arrayOfDouble2D getArrayOfDouble2D(int n, int m, double val)
{
    arrayOfDouble2D arrayDouble2D(n, vector<double>(m, val));
    return arrayDouble2D;
}

arrayOfDouble3D getArrayOfDouble3D(int n, int m, int p, double val)
{
    arrayOfDouble3D arrayDouble3D(n, vector<vector<double>>(m, vector<double>(p, val)));
    return arrayDouble3D;
}

arrayOfDouble4D getArrayOfDouble4D(int n, int m, int p, int q, double val)
{
    arrayOfDouble4D arrayDouble4D(n, vector<vector<vector<double>>>(m, vector<vector<double>>(p, vector<double>(q, val))));
    return arrayDouble4D;
}

double getMinimumOf3DArray(arrayOfDouble3D array3D, int& c)
{
    c = 0;
    double minValue = DBL_MAX;
    for (size_t d1 = 0; d1 < array3D.size(); d1++)
        for (size_t d2 = 0; d2 < array3D[d1].size(); d2++)
            for (size_t d3 = 0; d3 < array3D[d1][d2].size(); d3++)
            {
                minValue = min(minValue, array3D[d1][d2][d3]);
                if(minValue < array3D[d1][d2][d3] && minValue < 0.0)
                    c++;
            }
    return minValue;
}

double getMinimumOf2DArray(arrayOfDouble2D array2D)
{
    double minValue = DBL_MAX;
    for (size_t d1 = 0; d1 < array2D.size(); d1++)
        for (size_t d2 = 0; d2 < array2D[d1].size(); d2++)
        {
            minValue = min(minValue, array2D[d1][d2]);
        }
    return minValue;
}

double getMaximumOfArray(std::vector<double> vec)
{
    double maxValue = -DBL_MAX;
    for (auto v : vec)
        maxValue = max(maxValue, v);

    return maxValue;
}

double getMaximumOf2DArray(arrayOfDouble2D array2D)
{
    double maxValue = -DBL_MAX;
    for (size_t d1 = 0; d1 < array2D.size(); d1++)
        for (size_t d2 = 0; d2 < array2D[d1].size(); d2++)
        {
            maxValue = max(maxValue, array2D[d1][d2]);
        }
    return maxValue;
}

int getCountOfNegativeIn3DArray(arrayOfDouble3D array3D)
{
    int count = 0;
    for (size_t d1 = 0; d1 < array3D.size(); d1++)
        for (size_t d2 = 0; d2 < array3D[d1].size(); d2++)
            for (size_t d3 = 0; d3 < array3D[d1][d2].size(); d3++)
            {
                if (array3D[d1][d2][d3] < 0.0)
                    count++;
            }
    return count;
}

vector<double> linearize3DVector(arrayOfDouble3D array3D)
{
    // vector<double> data;
    size_t dim1 = array3D.size();
    size_t dim2 = array3D[0].size();
    size_t dim3 = array3D[0][0].size();
    vector<double> data(dim1 * dim2 * dim3, 0.0);
    for (size_t d1 = 0; d1 < dim1; d1++)
        for (size_t d2 = 0; d2 < dim2; d2++)
            for (size_t d3 = 0; d3 < dim3; d3++)
                data[d1 * dim2 * dim3 + d2 * dim3 + d3] = array3D[d1][d2][d3];

    return data;
}

string moreSigs(double d, int prec)
{
    std::stringstream ss;
    ss << d;
    ss.str("");
    ss << std::setprecision(prec) << std::fixed << d;

    std::string str;
    ss >> str;
    std::string::size_type s;
    for (s = str.length() - 1; s > 0; --s)
    {
        if (str[s] == '0')
            str.erase(s, 1);
        else
            break;
    }
    if (str[s] == '.')
        str.erase(s, 1);
    return str;
}

vector<string> listFiles(string path, string ext)
{
    bool gIgnoreHidden = true;
    string dotExt = "." + ext;
    vector<string> filelist;
    DIR *dirFile = opendir(path.c_str());
    if (dirFile)
    {
        struct dirent *hFile;
        errno = 0;
        while ((hFile = readdir(dirFile)) != NULL)
        {
            if (!strcmp(hFile->d_name, "."))
                continue;
            if (!strcmp(hFile->d_name, ".."))
                continue;

            // in linux hidden files all start with '.'
            if (gIgnoreHidden && (hFile->d_name[0] == '.'))
                continue;

            // dirFile.name is the name of the file. Do whatever string comparison
            // you want here. Something like:
            if (strstr(hFile->d_name, dotExt.c_str()))
                filelist.push_back(hFile->d_name);
            //cout<< "found an " << ext << " file " << hFile->d_name << endl;
        }
        closedir(dirFile);
    }
    return filelist;
}
