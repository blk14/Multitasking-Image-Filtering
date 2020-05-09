#include <iostream>
#include <fstream>
#include <mpi.h>
#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include <stdint.h>
#include <ctype.h>

using namespace std;

int INITIATOR = 0;
int description_size = 3 + 2 * sizeof(int) + 2 + sizeof(int) + 1 + 45;
int MAX_SEND_SIZE = 4000;


float Smooth[3][3] = {
    (float) 1/9, (float) 1/9, (float) 1/9,
    (float) 1/9, (float) 1/9, (float) 1/9,
    (float) 1/9, (float) 1/9, (float) 1/9
};

float Blur[3][3] = {
    (float)1/16, (float)2/16, (float)1/16,
    (float)2/16, (float)4/16, (float)2/16,
    (float)1/16, (float)2/16, (float)1/16
};

float Sharpen[3][3] = {
    0, (float)-2/3, 0,
    (float)-2/3, (float)11/3, (float)-2/3,
    0, (float)-2/3, 0
};

float MeanRemoval[3][3] = {
    -1, -1, -1,
    -1, 9, -1,
    -1, -1, -1
};

float Emboss[3][3] = {
    0, 1, 0,
    0, 0, 0,
    0, -1, 0
};

class Image {
 public:
    uint8_t** memblock;
    int width;
    int height;
    int maxVal;
    char type[3];

    void alloc() {
        this->memblock = new uint8_t*[this->height];
        for (int i = 0; i < this->height; i++) {
            this->memblock[i] = new uint8_t[this->width];
        }
    }
};

class Chunk {
 public:
    int height;
    int width;
    int maxVal;
    char type[3];
    uint8_t** memblock;

    void copy(Chunk* second) {
        second->height = height;
        second->width = width;
        second->maxVal = maxVal;
        second->alloc(this->height - 2, this->width - 6);
        for (int i = 0; i < this->height; i++) {
            for (int j = 0; j < this->width; j++) {
                second->memblock[i][j] = this->memblock[i][j];
            }
        }
    }

    void alloc(int size, int width) {         // + bordering
        this->height = size + 2;
        this->width = width + 6;
        this->memblock = new uint8_t*[this->height];
        for (int i = 0; i < this->height; i++) {
            this->memblock[i] = new uint8_t[this->width];
        }

        for (int i = 0; i < this->height; i++) {
            this->memblock[i][0] = 0;
            this->memblock[i][1] = 0;
            this->memblock[i][2] = 0;
            this->memblock[i][this->width - 1] = 0;
            this->memblock[i][this->width - 2] = 0;
            this->memblock[i][this->width - 3] = 0;
        }
        for (int j = 0; j < this->width; j++) {
            this->memblock[0][j] = 0;
            this->memblock[this->height - 1][j] = 0;
        }
    }
};

int minim(int a, int b) {
	if (a < b) {
		return a;
	}
	return b;
}

void splitSendRecv(bool ptype, uint8_t** memblock, int line_index, int width, int destination, int send_offset, int recv_offset) {
    int numberOfSR = width / MAX_SEND_SIZE;
    
    if (width % MAX_SEND_SIZE != 0) {
        numberOfSR += 1;
    }

    for (int i = 0; i < numberOfSR; i++) {
        if (i == numberOfSR - 1 && width % MAX_SEND_SIZE != 0) {
            if (ptype) {
                MPI_Send(&memblock[line_index][i * MAX_SEND_SIZE + send_offset], width - (i * MAX_SEND_SIZE), MPI_UINT8_T,
                        destination, 0, MPI_COMM_WORLD);
            } else {
                MPI_Recv(&memblock[line_index][i * MAX_SEND_SIZE + recv_offset], width - (i * MAX_SEND_SIZE), MPI_UINT8_T,
                        destination, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        } else {
            if (ptype) {
                MPI_Send(&memblock[line_index][i * MAX_SEND_SIZE + send_offset], MAX_SEND_SIZE, MPI_UINT8_T,
                        destination, 0, MPI_COMM_WORLD);
            } else {
                MPI_Recv(&memblock[line_index][i * MAX_SEND_SIZE + recv_offset], MAX_SEND_SIZE, MPI_UINT8_T,
                        destination, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
}
bool checkn(char buff[4]) {
    if (isdigit(buff[3])) {
        return true;
    }
    return false;
}

void readImage(char* argv[], Image* image, int nProcesses, MPI_Request request, MPI_Status status) {
    ifstream file (argv[1], ios::in|ios::binary);
    bool ok = true;
    if (file.is_open()) {
        file.read (image->type, 3);
        image->type[2] = '\0';
        char buff[4];
        file.read(buff, 1);
        while (buff[0] != '\n') {
            file.read(buff, 1);
        }

        file.read(buff, 4 * sizeof(char));
        ok = checkn(buff);
        image->width = atoi(buff);
        
        if (ok) {
            file.read(buff, 1);
            ok = true;
        }        
        file.read(buff, 4 * sizeof(char));
        ok = checkn(buff);
        image->height = atoi(buff);
        
        if (ok) {
            file.read(buff, 1);
            ok = true;
        }
        file.read(buff, sizeof(int));
        image->maxVal = atoi(buff);

    } else {
        cout << "maybe next time";
    }

    if (strcmp(image->type, "P6") == 0) {
        image->width *= 3;                     // AM SCHIMBAT DIMENSIUNEA PT P6!!!!
    }
    image->alloc();

    char buff[1];
    int step;
    if (nProcesses > 1) {
        step = image->height / (nProcesses - 1);
    } else {
        step = image->height;
    }
    int pCounter = 1;
    int previousStep = 0;
    int nextStep = step;
    int request_complete = 1;
    for (int i = 0; i < image->height; i++) {
        for (int j = 0; j < image->width; j++) {
            file.read(buff, 1);
            image->memblock[i][j] = buff[0];
        }
        if (nProcesses > 1) {
            if (strcmp(image->type, "P6") == 0) {
               splitSendRecv(true, image->memblock, i, image->width, pCounter, 0, 0);

            } else {
                MPI_Send(&image->memblock[i][0], image->width, MPI_UINT8_T, pCounter, 0, MPI_COMM_WORLD);
            }
            if (i == nextStep - 1) {
                pCounter ++;
                previousStep = nextStep;
                nextStep += step;
                if (nextStep < image->height && pCounter == nProcesses - 1) {
                    nextStep = image->height;
                }
            }
        }
    }
    file.close();
}


void printImage(char* argv[], Image* image) {
    ofstream outfile (argv[2], ios::out|ios::binary);
    if (outfile.is_open()) {
        outfile << image->type;
        outfile << "\n";
        if (strcmp(image->type, "P6") == 0) {
            image->width /= 3;
        }
        outfile << image->width << " " << image->height << "\n";
        if (strcmp(image->type, "P6") == 0) {
            image->width *= 3;
        }  
        outfile << image->maxVal << "\n";
    }
    for (int i = 0; i < image->height; i++) {
        for (int j = 0; j < image->width; j++) {
            outfile << image->memblock[i][j];
        }
    }
    outfile.close();
}


void getMyChunk(int rank, int nProcesses, Chunk* myChunk, int type, int height, int width, MPI_Request request, MPI_Status status) {
    int size = height / (nProcesses - 1);
    int request_complete = 0;

    if (rank == nProcesses - 1) {
        size = height - ((rank - 1) * size);
    }
    myChunk->alloc(size, width);

    for (int i = 0; i < size; i++) {
        if (type == 6) {
            splitSendRecv(false, myChunk->memblock, i + 1, width, 0, 0, 3);
        } else {
            MPI_Recv(&myChunk->memblock[i + 1][3], width, MPI_UINT8_T, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
}


vector<vector<float> > rotateMatrix(float buff[3][3]) {
    vector<vector<float> > kernel(3, vector<float>(3, 0));
    kernel[2][0] = buff[0][2];
    kernel[2][1] = buff[0][1];
    kernel[2][2] = buff[0][0];
    kernel[1][0] = buff[1][2];
    kernel[1][1] = buff[1][1];
    kernel[1][2] = buff[1][0];
    kernel[0][0] = buff[2][2];
    kernel[0][1] = buff[2][1];
    kernel[0][2] = buff[2][0];

    return kernel;
}


void applyFilters(Chunk* myChunk, int argc, char* argv[], int maxVal, int rank, int nProcesses, int type) {
    bool isBordered = false;  
    vector<vector<float> > bufferKernel(3, vector<float>(3,0));
    for (int ar = 3; ar < argc; ar++) {
         if (strcmp(argv[ar], "smooth") == 0) {
            bufferKernel = rotateMatrix(Smooth);
        } else if (strcmp(argv[ar], "blur") == 0) {
            bufferKernel = rotateMatrix(Blur);
        } else if (strcmp(argv[ar], "sharpen") == 0) {
            bufferKernel = rotateMatrix(Sharpen);
        } else if (strcmp(argv[ar], "mean") == 0) {
            bufferKernel = rotateMatrix(MeanRemoval);
        } else if (strcmp(argv[ar], "emboss") == 0) {
            bufferKernel = rotateMatrix(Emboss);
        } else {
            cout << "Something strange with your args\n";
        }
        if (rank == 1 && nProcesses > 2) {
            if (type == 6) {
                
                splitSendRecv(true, myChunk->memblock, myChunk->height - 2, myChunk->width, rank + 1, 0, 0);

                splitSendRecv(false, myChunk->memblock, myChunk->height - 1, myChunk->width, rank + 1, 0, 0);

            } else {
                MPI_Send(&myChunk->memblock[myChunk->height - 2][3], myChunk->width - 6,
                         MPI_UINT8_T, rank + 1, 0, MPI_COMM_WORLD);

                MPI_Recv(&myChunk->memblock[myChunk->height - 1][3], myChunk->width - 6, MPI_UINT8_T,
                        rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        } else if (rank == nProcesses - 1 && nProcesses > 2) {
            if (type == 6) {

                splitSendRecv(true, myChunk->memblock, 1, myChunk->width, rank - 1, 0, 0);

                splitSendRecv(false, myChunk->memblock, 0, myChunk->width, rank - 1, 0, 0);

            } else {
                MPI_Send(&myChunk->memblock[1][3], myChunk->width - 6,
                         MPI_UINT8_T, rank - 1, 0, MPI_COMM_WORLD);

                MPI_Recv(&myChunk->memblock[0][3], myChunk->width - 6, MPI_UINT8_T,
                        rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        } else if (rank != 1 && rank != nProcesses - 1) {
            if (type == 6) {

                splitSendRecv(true, myChunk->memblock, myChunk->height - 2, myChunk->width, rank + 1, 0, 0);
                splitSendRecv(true, myChunk->memblock, 1, myChunk->width, rank - 1, 0, 0);
                
                splitSendRecv(false, myChunk->memblock, myChunk->height - 1, myChunk->width, rank + 1, 0, 0);
                splitSendRecv(false, myChunk->memblock, 0, myChunk->width, rank - 1, 0, 0);

            } else {
                MPI_Send(&myChunk->memblock[myChunk->height - 2][3], myChunk->width - 6,
                         MPI_UINT8_T, rank + 1, 0, MPI_COMM_WORLD);

                MPI_Send(&myChunk->memblock[1][3], myChunk->width - 6,
                         MPI_UINT8_T, rank - 1, 0, MPI_COMM_WORLD);

                MPI_Recv(&myChunk->memblock[myChunk->height - 1][3], myChunk->width - 6, MPI_UINT8_T,
                        rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                MPI_Recv(&myChunk->memblock[0][3], myChunk->width - 6, MPI_UINT8_T,
                        rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        Chunk* initialChunk = new Chunk();
        myChunk->copy(initialChunk);

        for (int i = 1; i < myChunk->height - 1; i++) {
            for (int j = 3; j < myChunk->width - 3; j++) {
                float sum = 0;
                if (type == 6) {
                    sum += (float)initialChunk->memblock[i - 1][j - 3] * bufferKernel[0][0] +
                            (float)initialChunk->memblock[i - 1][j] * bufferKernel[0][1] +
                            (float)initialChunk->memblock[i - 1][j + 3] * bufferKernel[0][2] +
                            (float)initialChunk->memblock[i][j - 3] * bufferKernel[1][0] +
                            (float)initialChunk->memblock[i][j] * bufferKernel[1][1] +
                            (float)initialChunk->memblock[i][j + 3] * bufferKernel[1][2] +
                            (float)initialChunk->memblock[i + 1][j - 3] * bufferKernel[2][0] +
                            (float)initialChunk->memblock[i + 1][j] * bufferKernel[2][1] +
                            (float)initialChunk->memblock[i + 1][j + 3] * bufferKernel[2][2];
                } else {
                    sum += (float)initialChunk->memblock[i - 1][j - 1] * bufferKernel[0][0] +
                            (float)initialChunk->memblock[i - 1][j] * bufferKernel[0][1] +
                            (float)initialChunk->memblock[i - 1][j + 1] * bufferKernel[0][2] +
                            (float)initialChunk->memblock[i][j - 1] * bufferKernel[1][0] +
                            (float)initialChunk->memblock[i][j] * bufferKernel[1][1] +
                            (float)initialChunk->memblock[i][j + 1] * bufferKernel[1][2] +
                            (float)initialChunk->memblock[i + 1][j - 1] * bufferKernel[2][0] +
                            (float)initialChunk->memblock[i + 1][j] * bufferKernel[2][1] +
                            (float)initialChunk->memblock[i + 1][j + 1] * bufferKernel[2][2];
                }
                if (sum > maxVal) {
                    sum = maxVal;
                }
                if (sum < 0) {
                    sum = 0;
                }
                myChunk->memblock[i][j] = (uint8_t) sum;
            }
        }
        for (int i = 0; i < initialChunk->height; i++) {
            delete[] initialChunk->memblock[i];
        }
        delete[] initialChunk->memblock;
        delete initialChunk;
    }
}

void recollectImage(Image* image, MPI_Status status, int height, int width, int nProcesses, int type) {
    int line_index; 
    for (int i = 0; i < height; i++) {
        if (type == 6) {
            MPI_Recv(&line_index, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

            splitSendRecv(false, image->memblock, line_index, image->width, status.MPI_SOURCE, 0, 0);
        } else {
            MPI_Recv(&line_index, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

            MPI_Recv(&image->memblock[line_index][0], width, MPI_UINT8_T,
                    status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
}


void sendChunk(Chunk* myChunk, int type, int rank, int nProcesses, int height, int width) {
    int offset = (rank - 1) * (height / (nProcesses - 1));
    int line_index;
    for (int i = 1; i < myChunk->height - 1; i++) {
        if (type == 6) {
            line_index = offset + i - 1;
            MPI_Send(&line_index, 1, MPI_INT, INITIATOR, 0, MPI_COMM_WORLD);

            splitSendRecv(true, myChunk->memblock, i, width, INITIATOR, 3, 0);
        } else {
            line_index = offset + i - 1;
            MPI_Send(&line_index, 1, MPI_INT, INITIATOR, 0, MPI_COMM_WORLD);

            MPI_Send(&myChunk->memblock[i][3], width, MPI_UINT8_T,
                    INITIATOR, 0, MPI_COMM_WORLD);
        }
    }
}

int main(int argc, char* argv[]) {
    int rank;
	int nProcesses;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);
    MPI_Request request;
    MPI_Status  status;
    int type = 0;                    // 5 or 6
    int height = 0;
    int width = 0;
    int maxVal = 0;

    Image* image = new Image();
    if (rank == INITIATOR) {
        readImage(argv, image, nProcesses, request, status);
        if (strcmp(image->type, "P6") == 0) {
            type = 6;
        } else {
            type = 5;
        }
        height = image->height;
        width = image->width;
        maxVal = image->maxVal;
        if (nProcesses == 1) {
            Chunk* myChunk = new Chunk();
            myChunk->alloc(height, width);
            for (int i = 0; i < height; i++) {
                for (int j = 0; j < width; j++) {
                    myChunk->memblock[i + 1][j + 3] = image->memblock[i][j];
                }
            }
            applyFilters(myChunk, argc, argv, maxVal, rank, nProcesses, type);
            for (int i = 0; i < height; i++) {
                for (int j = 0; j < width; j++) {
                    image->memblock[i][j] = myChunk->memblock[i + 1][j + 3];
                }
            }
            printImage(argv, image);

            for (int i = 0; i < image->height; i++) {
                delete[] image->memblock[i];
            }
            delete[] image->memblock;
            delete image;

            for (int i = 0; i < myChunk->height; i++) {
                delete[] myChunk->memblock[i];
            }
            delete[] myChunk->memblock;
            delete myChunk;
            
            MPI_Finalize();
            return 0;
        }
    }
    
    MPI_Bcast(&type, 1, MPI_INT, INITIATOR, MPI_COMM_WORLD);    // find out the type of the image
    MPI_Bcast(&height, 1, MPI_INT, INITIATOR, MPI_COMM_WORLD);    // find out the height of the image
    MPI_Bcast(&width, 1, MPI_INT, INITIATOR, MPI_COMM_WORLD);    // find out the width of the image (included P6)
    MPI_Bcast(&maxVal, 1, MPI_INT, INITIATOR, MPI_COMM_WORLD);    // find out the maxVal of the image 

    
    Chunk* myChunk = new Chunk();
    if (rank != INITIATOR) {
        getMyChunk(rank, nProcesses, myChunk, type, height, width, request, status);
        applyFilters(myChunk, argc, argv, maxVal, rank, nProcesses, type);
    }

    if (rank == INITIATOR) {
        recollectImage(image, status, height, width, nProcesses, type);
        printImage(argv, image);

        for (int i = 0; i < image->height; i++) {
            delete[] image->memblock[i];
        }
        delete[] image->memblock;
        delete image;

    } else {
        sendChunk(myChunk, type, rank, nProcesses, height, width);

        for (int i = 0; i < myChunk->height; i++) {
            delete[] myChunk->memblock[i];
        }
        delete[] myChunk->memblock;
        delete myChunk;
    }

    MPI_Finalize();
    return 0;
}