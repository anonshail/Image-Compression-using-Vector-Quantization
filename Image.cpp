//*****************************************************************************
//
// Image.cpp : Defines the class operations on images
//
// Author - Parag Havaldar
// Code used by students as starter code to display and modify images
//
//*****************************************************************************

#include "Image.h"
#include <vector>

using namespace std;

// Constructor and Desctructors
MyImage::MyImage() 
{
	Data = NULL;
	OneChannelData = NULL;
	Width = -1;
	Height = -1;
	ImagePath[0] = 0;
}

MyImage::~MyImage()
{
	if ( Data )
		delete Data;
	if ( OneChannelData )
		delete OneChannelData;
}


// Copy constructor
MyImage::MyImage( MyImage *otherImage)
{
	Height = otherImage->Height;
	Width  = otherImage->Width;
	Data   = new char[Width*Height*3];
	strcpy(otherImage->ImagePath, ImagePath );

	for ( int i=0; i<(Height*Width*3); i++ )
	{
		Data[i]	= otherImage->Data[i];
	}

	if (otherImage->OneChannelData) {
		OneChannelData = new char[Width * Height];
		for (int i = 0; i < (Height * Width); i++)
		{
			OneChannelData[i] = otherImage->OneChannelData[i];
		}
	}


}



// = operator overload
MyImage & MyImage::operator= (const MyImage &otherImage)
{
	Height = otherImage.Height;
	Width  = otherImage.Width;
	Data   = new char[Width*Height*3];
	strcpy( (char *)otherImage.ImagePath, ImagePath );

	for ( int i=0; i<(Height*Width*3); i++ )
	{
		Data[i]	= otherImage.Data[i];
	}

	if (otherImage.OneChannelData) {
		OneChannelData = new char[Width * Height];
		for (int i = 0; i < (Height * Width); i++)
		{
			OneChannelData[i] = otherImage.OneChannelData[i];
		}
	}
	
	return *this;

}


// MyImage::ReadImage
// Function to read the image given a path
bool MyImage::ReadImage()
{

	// Verify ImagePath
	if (ImagePath[0] == 0 || Width < 0 || Height < 0 )
	{
		fprintf(stderr, "Image or Image properties not defined");
		fprintf(stderr, "Usage is `Image.exe Imagefile w h`");
		return false;
	}
	
	// Create a valid output file pointer
	FILE *IN_FILE;
	IN_FILE = fopen(ImagePath, "rb");
	if ( IN_FILE == NULL ) 
	{
		fprintf(stderr, "Error Opening File for Reading");
		return false;
	}

	// Create and populate RGB buffers
	int i;
	char *Rbuf = new char[Height*Width]; 
	char *Gbuf = new char[Height*Width]; 
	char *Bbuf = new char[Height*Width]; 

	for (i = 0; i < Width*Height; i ++)
	{
		Rbuf[i] = fgetc(IN_FILE);
	}
	for (i = 0; i < Width*Height; i ++)
	{
		Gbuf[i] = fgetc(IN_FILE);
	}
	for (i = 0; i < Width*Height; i ++)
	{
		Bbuf[i] = fgetc(IN_FILE);
	}
	
	// Allocate Data structure and copy
	Data = new char[Width*Height*3];
	for (i = 0; i < Height*Width; i++)
	{
		Data[3*i]	= Bbuf[i];
		Data[3*i+1]	= Gbuf[i];
		Data[3*i+2]	= Rbuf[i];
	}

	// Clean up and return
	delete Rbuf;
	delete Gbuf;
	delete Bbuf;
	fclose(IN_FILE);

	return true;

}

// MyImage::ReadImageOneChannel
// Function to read the image given a path in one channel mode
bool MyImage::ReadImageOneChannel()
{

	// Verify ImagePath
	if (ImagePath[0] == 0 || Width < 0 || Height < 0)
	{
		fprintf(stderr, "Image or Image properties not defined");
		fprintf(stderr, "Usage is `Image.exe Imagefile w h`");
		return false;
	}

	// Create a valid output file pointer
	FILE* IN_FILE;
	IN_FILE = fopen(ImagePath, "rb");
	if (IN_FILE == NULL)
	{
		fprintf(stderr, "Error Opening File for Reading");
		return false;
	}

	// Create and populate RGB buffers
	int i;
	char* Rbuf = new char[Height * Width];
	char* Gbuf = new char[Height * Width];
	char* Bbuf = new char[Height * Width];

	for (i = 0; i < Width * Height; i++)
	{
		Rbuf[i] = fgetc(IN_FILE);
		Gbuf[i] = Rbuf[i];
		Bbuf[i] = Rbuf[i];
	}

	// Allocate Data structure and copy
	Data = new char[Width * Height * 3];
	OneChannelData = new char[Width * Height];
	for (i = 0; i < Height * Width; i++)
	{
		Data[3 * i] = Bbuf[i];
		OneChannelData[i] = Bbuf[i];
		Data[3 * i + 1] = Gbuf[i];
		Data[3 * i + 2] = Rbuf[i];
	}

	// Clean up and return
	delete Rbuf;
	delete Gbuf;
	delete Bbuf;
	fclose(IN_FILE);

	return true;

}



// MyImage functions defined here
bool MyImage::WriteImage()
{
	// Verify ImagePath
	// Verify ImagePath
	if (ImagePath[0] == 0 || Width < 0 || Height < 0 )
	{
		fprintf(stderr, "Image or Image properties not defined");
		return false;
	}
	
	// Create a valid output file pointer
	FILE *OUT_FILE;
	OUT_FILE = fopen(ImagePath, "wb");
	if ( OUT_FILE == NULL ) 
	{
		fprintf(stderr, "Error Opening File for Writing");
		return false;
	}

	// Create and populate RGB buffers
	int i;
	char *Rbuf = new char[Height*Width]; 
	char *Gbuf = new char[Height*Width]; 
	char *Bbuf = new char[Height*Width]; 

	for (i = 0; i < Height*Width; i++)
	{
		Bbuf[i] = Data[3*i];
		Gbuf[i] = Data[3*i+1];
		Rbuf[i] = Data[3*i+2];
	}

	
	// Write data to file
	for (i = 0; i < Width*Height; i ++)
	{
		fputc(Rbuf[i], OUT_FILE);
	}
	for (i = 0; i < Width*Height; i ++)
	{
		fputc(Gbuf[i], OUT_FILE);
	}
	for (i = 0; i < Width*Height; i ++)
	{
		fputc(Bbuf[i], OUT_FILE);
	}
	
	// Clean up and return
	delete Rbuf;
	delete Gbuf;
	delete Bbuf;
	fclose(OUT_FILE);

	return true;

}



// =================== UTILITY FUNCTIONS FOR VECTOR QUANTIZATION =========================
double distanceFormula(const vector<int>& v1, const vector<int>& v2) {
	double s = 0;
	for (int i = 0; i < v1.size(); i++) {
		s += ((v2[i] - v1[i]) * (v2[i] - v1[i]));
	}
	double ans = sqrt(s);
	return ans;
}

void clusterIt(const vector<vector<int>>& codewords, const int sizeOfVector, vector<vector<vector<int>>>& clusters, int Width, int Height, const char* Data, vector<pair<bool, int>>& updateTable) {
	int W = int(ceil(sqrt(sizeOfVector))), H = sizeOfVector / W, minUp = int(clusters.size() / 13);
	
	if (minUp < 10) {
		minUp = 10;
		if (sqrt(sizeOfVector) >= 2) {
			minUp *= (W/2);
		}
		if (minUp > 20) minUp = 20;
	}
	printf("minUp: %d\n", minUp);
	if (clusters[0].size() <= 0) {
		for (int i = 0; i < Height / W; i++) {
			for (int j = 0; j < Width / H; j++) {
				vector<int> v;
				v.reserve(sizeOfVector);

				for (int row = H * i; H * (1 + i) > row; row++) {
					for (int col = W * j; col < W * (j + 1); col++) {
						v.push_back(Data[col + Width * row]);
					}
				}

				int minIndex = 0;
				double min = distanceFormula(codewords[minIndex], v);

				for (int i = 0; i < codewords.size(); i++) {
					double d = distanceFormula(codewords[i], v);
					if (d < min) {
						minIndex = i;
						min = d;
					}
				}

				clusters[minIndex].push_back(v);
			}
		}
	}
	else {
		int sMin = INT_MAX;

		// setting all clusters to not update, change when updated
		for (int i = 0; i < clusters.size(); i++) {
			updateTable[i].first = false;
		}

		// local clusters
		vector<vector<vector<int>>> lClusters(clusters.size());

		for (int i = 0; i < clusters.size(); i++) {
			if (updateTable[i].second >= minUp) {
				continue;
			}
			if (i < sMin) {
				sMin = i;
			}

			for (int j = 0; clusters[i].size() > j; j++) {
				vector<int>& v = clusters[i][j];
				double md = distanceFormula(codewords[sMin], v);
				int minIndex = sMin;
				for (int k = minIndex + 1; k < codewords.size(); k++) {
					double d = distanceFormula(codewords[k], v);
					if (updateTable[k].second >= minUp) {
						continue;
					} else if (d < md) {
						minIndex = k;
						md = d;
					}
				}
				lClusters[minIndex].push_back(v);

				if (!updateTable[minIndex].first && lClusters[minIndex].size() - 1 < clusters[minIndex].size() && clusters[minIndex][lClusters[minIndex].size() - 1] == v) {
					continue;
				} else {
					updateTable[minIndex].second = 0;
					updateTable[minIndex].first = true;
				}
			}
		}

		for (int i = 0; clusters.size() > i; i++) {
			if (!updateTable[i].first) {
				updateTable[i].second++;
			} else {
				clusters[i] = lClusters[i];
			}
		}
	}
}

int updatedCodewords(const vector<vector<vector<int>>>& clusters, const int sizeOfVector, vector<vector<int>>& codewords, vector<pair<bool, int>>& updateTable) {
	int updates = 0;

	for (int i = 0; i < clusters.size(); i++) {
		if (updateTable[i].first != false) {
			vector<int> v(sizeOfVector);

			if (clusters[i].size() == 0) {
				// Incase clusters are empty, assign random codewords
				srand(time(NULL));
				for (int j = 0; j < sizeOfVector; j++) {
					v[j] = rand() % 256;
				}
				updates++;
				codewords[i] = v;
			} else {
				for (int j = 0; j < sizeOfVector; j++) {
					double s = 0.0;
					for (auto vtr : clusters[i]) {
						s = s + (1 * vtr[j] / double(clusters[i].size()));
					}
					v[j] = int(s);
				}
				updates++;
				codewords[i] = v;
			}
		}
	}

	return updates;
}

void quantizeData(const vector<vector<int>>& codewords, const int sizeOfVector, const char* Data, const int Width, const int Height, vector<vector<int>>& quantData) {
	int W = int(ceil(sqrt(sizeOfVector))), H = sizeOfVector / W, qW = int(Width / W), qH = int(Height / H);
	quantData.reserve(qH);

	printf("W = %d\nH = %d\nWidth = %d\nHeight = %d\nqW = %d\nqH = %d\n", W, H, Width, Height, qW, qH);

	for (int i = 0; i < qH; i++) {
		vector<int> row(qW);

		for (int j = 0; j < qW; j++) {
			int minIndex = 0;
			vector<int> v;
			v.reserve(sizeOfVector);

			for (int row = H * i; row < H * (i + 1); row++) {
				for (int col = W * j; col < W * (j + 1); col++) {
					v.push_back(Data[col + Width * row]);
				}
			}

			double min = distanceFormula(codewords[minIndex], v);
			for (int i = 1; i < codewords.size(); i++) {
				double d = distanceFormula(codewords[i], v);
				if (d < min) {
					min = d;
					minIndex = i;
				}
			}

			row[j] = minIndex;
		}

		quantData.push_back(row);
	}
}

void VQEncode(const int numOfClusters, const int sizeOfVector, vector<vector<int>>& codewords, vector<vector<vector<int>>>& clusters, vector<vector<int>>& quantData, vector<pair<bool, int>>& updateTable, int Width, int Height, const char* Data) {
	int numUpdates = numOfClusters, count = 0, maxClusterIterations = 1000;

	while (count < maxClusterIterations && numUpdates > 1) {
		clusterIt(codewords, sizeOfVector, clusters, Width, Height, Data, updateTable);
		numUpdates = updatedCodewords(clusters, sizeOfVector, codewords, updateTable);

		printf("Iteration %d: No. of clusters updated: %d\n", ++count, numUpdates);
	}

	quantizeData(codewords, sizeOfVector, Data, Width, Height, quantData);
}

void VQDecode(const vector<vector<int>>& codewords, const vector<vector<int>>& quantData, const int sizeOfVector, const int Width, const int Height, char* Data) {
	int W = int(ceil(sqrt(sizeOfVector))), H = sizeOfVector / W;

	for (int i = 0; i < quantData.size(); i++) {
		for (int j = 0; j < quantData[i].size(); j++) {
			int indexer = quantData[i][j];
			vector<int> v = codewords[indexer];
			vector<int>::iterator it = v.begin();

			for (int row = H * i; row < H * (i + 1); row++) {
				for (int col = W * j; col < W * (j + 1); col++) {
					Data[row * Width + col] = byte(*it);
					it++;
				}
			}

		}
	}
}

// This function compresses image data using Vector Quantization
bool MyImage::Modify(const int numOfClusters, const int sizeOfVector)
{	
	if (!OneChannelData) return false;

	// Create VQ Resources and Reserve Memory
	vector<vector<int>> codewords;				// codebook vectors
	vector<vector<vector<int>>> clusters;		// cluster data
	vector<vector<int>> quantData;				// result of quantized data after VQ
	vector<pair<bool, int>> updateTable;		// Update Table which indicates if the cluster is changed

	// Allocate memory
	codewords.reserve(numOfClusters);	 
	clusters.reserve(numOfClusters);
	updateTable.reserve(numOfClusters);

	// Initialize first N vectors by dividing the diagnal equally
	for (int i = 0; i < numOfClusters; i++) {
		vector<int> v;
		v.reserve(sizeOfVector);

		for (int j = 0; j < sizeOfVector; j++) {
			v.push_back(int(256.0 * (i + 1) / double(numOfClusters) + 1));
		}

		// Add empty cluster to clusters
		clusters.push_back(vector<vector<int>>());
		// Add code word for that cluster
		codewords.push_back(v);
		// Updataing
		updateTable.push_back(pair<bool, int>(true, 0));
	}


	// Encode image using Vector Quantization 
	VQEncode(numOfClusters, sizeOfVector, codewords, clusters, quantData, updateTable, Width, Height, OneChannelData);

	// Decode image and set the data as required :)
	char* compressedOneChannel = new char[Width * Height];
	VQDecode(codewords, quantData, sizeOfVector, Width, Height, compressedOneChannel);

	char* compressedThreeChannelData = new char[Width * Height * 3];
	for (int i = 0; i < Height * Width; i++)
	{
		compressedThreeChannelData[3 * i] = compressedOneChannel[i];
		compressedThreeChannelData[3 * i + 1] = compressedOneChannel[i];
		compressedThreeChannelData[3 * i + 2] = compressedOneChannel[i];
	}

	char* oldData = getOneChannelImageData();
	delete oldData;
	setOneChannelImageData(compressedOneChannel);
	setImageData(compressedThreeChannelData);

	return false;
}