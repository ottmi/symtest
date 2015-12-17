#include "heatmap.h"

//==============================================================//
// Color-coding related settings and functions
//==============================================================//

static int HEATMAP_COLOR_NUM = 11;
static string HEATMAP_COLOR_MAP[] = {"#AA0000","#FF2A2A", "#FF6600", "#FF9955", "#FFCC00", "#DDFF55", "#99FF55", "#66FF00", "#00AA00", "#005500", "#002B00"};
static string HEATMAP_COLOR_DESC[] = {"= 0.0", "&lt; 0.1", "&lt; 0.2", "&lt; 0.3", "&lt; 0.4", "&lt; 0.5", "&lt; 0.6", "&lt; 0.7", "&lt; 0.8", "&lt; 0.9", "&#x2264; 1.0"};

// get the corresponding color code according to the value
string getColor(double value) {
	
	if (value == 0.0)
		return HEATMAP_COLOR_MAP[0];
	else
		return HEATMAP_COLOR_MAP[(int)(value*10.0) + 1];
}

//==============================================================//
// Other functions for SVG generation
//==============================================================//

#define HEATMAP_BLOCK_DIM 5 // the dimension of the squares for heatmap
#define HEATMAP_BLOCK_DIM_FEW_SEQ 10 // the dimension of the squares for heatmap if the number of sequences is "FEW"
#define HEATMAP_FEW_SEQ_NUM 10 // the number of sequences is regarded as "FEW" if the number of sequences <= 10
#define HEATMAP_FONT "Arial" // Before: "Bitstream Vera Sans"
#define HEATMAP_FONT_SIZE 4 // the font size of the text for heatmap
#define HEATMAP_FONT_SIZE_FEW_SEQ 8 // the font size of the text for heatmap if the number of sequences is "FEW"
#define HEATMAP_TEXT_ALIGN "end" // the alignment of the text for heatmap
#define HEATMAP_LEGEND_BLOCK_DIM 10 // the dimension of the squares for the legend of the heatmap
#define HEATMAP_LEGEND_FONT "Arial" // Before: "Bitstream Vera Sans"
#define HEATMAP_LEGEND_FONT_SIZE 8 // the font size of the text for the legend of the heatmap
#define HEATMAP_LEGEND_TEXT_ALIGN "front" // the alignment of the text for the legend of the heatmap

#define HEATMAP_GRAY_COLOR "#cccccc" // gray color code for heatmap
#define HEATMAP_NAME_ANGLE 315 // angle of the sequence name for heatmap
#define HEATMAP_LEGEND_DESC_ANGLE 0 // angle of the description for the legend of the heatmap

#define HEATMAP_TOP_BORDER 5
#define HEATMAP_BOTTOM_BORDER 5
#define HEATMAP_LEFT_BORDER 5
#define HEATMAP_RIGHT_BORDER 5
#define HEATMAP_GAP_BW_TEXT 3
#define HEATMAP_GAP_BW_LEGEND 20
#define HEATMAP_LEGEND_GAP_BW_DESC 5

// Output the SVG Header
void outputSVGHeader(ofstream& fout, int width, int height) {
	fout << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" << endl;
	fout << "<svg" << endl;
	fout << "\txmlns:svg=\"http://www.w3.org/2000/svg\"" << endl;
	fout << "\txmlns=\"http://www.w3.org/2000/svg\"" << endl;
	fout << "\tversion=\"1.0\"" << endl;
	fout << "\twidth=\"" << width << "\"" << endl;
	fout << "\theight=\"" << height << "\"" << endl;
	fout << "\tid=\"svg2\">" << endl;
	fout << "\t<defs id=\"defs4\" />" << endl;
}

// Output the SVG Tail
void outputSVGTail(ofstream& fout) {
	fout << "</svg>" << endl;
}
// Output the square
void outputSquare(ofstream& fout, double x, double y, string color, int size) {
	fout << "<rect" << endl;
	fout << "\twidth=\"" << size << "\"" << endl;
	fout << "\theight=\"" << size << "\"" << endl;
	fout << "\tx=\"" << x << "\"" << endl;
	fout << "\ty=\"" << y << "\"" << endl;
	fout << "\tstyle=\"fill:" << color << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\" />" << endl;
}

// output the sequence name
void outputSeqName(ofstream& fout, double x, double y, int angle, string txt, string font, int fontSize, string text_align) {
	fout << "<text" << endl;
	fout << "\tx=\"" << x << "\" y=\"" << y << "\"" << endl;
	fout << "\tstyle=\"font-size:" << fontSize << "px;font-style:normal;font-weight:normal;text-align:" << text_align << ";text-anchor:" << text_align << ";fill:black;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;font-family:" << font << "\"" << endl;
	if (angle > 0)
		fout << "\ttransform=\"rotate(" << angle << " " << x << " " << y << ")\"" << endl;
	fout << "\t>" << txt << "</text>" << endl;
}

// get the longest length among the first "num" strings inside the array
int longestLenFirst(string* seqArray, int num, int totNum) {
	if (totNum == 0)
		return 0;
	unsigned int max_len = seqArray[0].length();
	for (int i=1; i<num && i<totNum; i++)
		if (seqArray[i].length() > max_len)
			max_len = seqArray[i].length();
	return (int) max_len;
}

// get the longest length among the last "num" strings inside the array
int longestLenLast(string* seqArray, int num, int totNum) {
	if (totNum == 0)
		return 0;
	unsigned int max_len = seqArray[totNum-1].length();
	int i;
	for (i=totNum-2; i>=totNum-num && i>=0; i--)
		if (seqArray[i].length() > max_len)
			max_len = seqArray[i].length();
	return (int) max_len;
}

// get the longest length among all the strings inside the array
int longestLen(string* seqArray, int totNum) {
	if (totNum == 0)
		return 0;
	unsigned int max_len = seqArray[0].length();
	int i;
	for (i=1; i<totNum; i++)
		if (seqArray[i].length() > max_len)
			max_len = seqArray[i].length();
	return (int) max_len;
}

// output the legend
void outputLegend(ofstream& fout, int x, int y) {

	double x_text, y_text;
	x_text = x + HEATMAP_LEGEND_BLOCK_DIM + HEATMAP_LEGEND_GAP_BW_DESC;
	y_text = y + HEATMAP_LEGEND_BLOCK_DIM*0.8;
	int i;
	for (i=0; i<HEATMAP_COLOR_NUM; i++) {
		// print out the square
		outputSquare(fout, x, y, HEATMAP_COLOR_MAP[i], HEATMAP_LEGEND_BLOCK_DIM);
		// print out the description for each color
		outputSeqName(fout, x_text, y_text, HEATMAP_LEGEND_DESC_ANGLE, HEATMAP_COLOR_DESC[i], HEATMAP_LEGEND_FONT, HEATMAP_LEGEND_FONT_SIZE, HEATMAP_LEGEND_TEXT_ALIGN);
		y += HEATMAP_LEGEND_BLOCK_DIM;
		y_text += HEATMAP_LEGEND_BLOCK_DIM;
	}
}

// Output the triangular heatmap
void outputTriHeatmap(string prefixOut, string* seqNames, double* cij, int seqNum) {
	
	int block_size = HEATMAP_BLOCK_DIM;
	int font_size = HEATMAP_FONT_SIZE;
	if (seqNum <= HEATMAP_FEW_SEQ_NUM) {
		block_size = HEATMAP_BLOCK_DIM_FEW_SEQ;
		font_size = HEATMAP_FONT_SIZE_FEW_SEQ;
	}

	int longestFirstTextLen = longestLenFirst(seqNames, 3, seqNum);
	int longestLastTextLen = longestLenLast(seqNames, 3, seqNum);
	int width = HEATMAP_LEFT_BORDER + longestFirstTextLen*font_size/2.0 + HEATMAP_GAP_BW_TEXT + 
				block_size*(seqNum-1) + HEATMAP_GAP_BW_LEGEND + HEATMAP_LEGEND_BLOCK_DIM +
				HEATMAP_LEGEND_GAP_BW_DESC + HEATMAP_LEGEND_FONT_SIZE*2.5 + HEATMAP_RIGHT_BORDER;
	int height = HEATMAP_TOP_BORDER + block_size*seqNum + (longestLastTextLen*font_size/2.0 - block_size) +
				 HEATMAP_BOTTOM_BORDER;
	int height_legend = HEATMAP_TOP_BORDER + HEATMAP_LEGEND_BLOCK_DIM * HEATMAP_COLOR_NUM + HEATMAP_BOTTOM_BORDER;
	if (height < height_legend)
		height = height_legend;
	double y = HEATMAP_TOP_BORDER;
	double x;
	double y_coord, x_coord;
	int i, j;
	
	string outFileName = prefixOut;
	outFileName.append(".triangular.svg");
	ofstream fout;
	fout.open(outFileName.c_str());

	// output the header
	outputSVGHeader(fout, width, height);
	
	// output the triangular matrix
	for (i=0; i<seqNum; i++) {
		
		x = HEATMAP_LEFT_BORDER + longestFirstTextLen*font_size/2.0 + HEATMAP_GAP_BW_TEXT + i*block_size;
		
		// print out the sequence name
		if (seqNames[i].length() > 0) {
			x_coord = x - HEATMAP_GAP_BW_TEXT;
			y_coord = y + block_size*0.8;
			outputSeqName(fout, x_coord, y_coord, HEATMAP_NAME_ANGLE, seqNames[i], HEATMAP_FONT, font_size, HEATMAP_TEXT_ALIGN);
		}
		
		for (j=i+1; j<seqNum; j++) {
			// print out the squares representing the Cij values
			outputSquare(fout, x, y, getColor(cij[i*seqNum+j]), block_size);
			x += block_size;
		}
		
		y += block_size;
	}
	
	// output the legend
	y = HEATMAP_TOP_BORDER;
	x = HEATMAP_LEFT_BORDER + longestFirstTextLen*font_size/2.0 + HEATMAP_GAP_BW_TEXT + 
				block_size*(seqNum-1) + HEATMAP_GAP_BW_LEGEND;
	outputLegend(fout, x, y);
	
	// output the tail
	outputSVGTail(fout);

	fout.close();
}

// Output the full heatmap
void outputFullHeatmap(string prefixOut, string* seqNames, double* cij, int seqNum) {
	
	int block_size = HEATMAP_BLOCK_DIM;
	int font_size = HEATMAP_FONT_SIZE;
	if (seqNum <= HEATMAP_FEW_SEQ_NUM) {
		block_size = HEATMAP_BLOCK_DIM_FEW_SEQ;
		font_size = HEATMAP_FONT_SIZE_FEW_SEQ;
	}
	
	int longestTextLen = longestLen(seqNames, seqNum);
	int width = HEATMAP_LEFT_BORDER + longestTextLen*font_size/2.0 + HEATMAP_GAP_BW_TEXT + 
				block_size*seqNum + HEATMAP_GAP_BW_LEGEND + HEATMAP_LEGEND_BLOCK_DIM +
				HEATMAP_LEGEND_GAP_BW_DESC + HEATMAP_LEGEND_FONT_SIZE*2.5 + HEATMAP_RIGHT_BORDER;
	int height = HEATMAP_TOP_BORDER + block_size*seqNum + HEATMAP_BOTTOM_BORDER;
	int height_legend = HEATMAP_TOP_BORDER + HEATMAP_LEGEND_BLOCK_DIM * HEATMAP_COLOR_NUM + HEATMAP_BOTTOM_BORDER;
	if (height < height_legend)
		height = height_legend;
	double y = HEATMAP_TOP_BORDER;
	double x;
	double y_coord, x_coord;
	int i, j;
	
	string outFileName = prefixOut;
	outFileName.append(".full.svg");
	ofstream fout;
	fout.open(outFileName.c_str());

	// output the header
	outputSVGHeader(fout, width, height);
	
	// output the full matrix
	for (i=0; i<seqNum; i++) {
		
		x = HEATMAP_LEFT_BORDER + longestTextLen*font_size/2.0 + HEATMAP_GAP_BW_TEXT;
		
		// print out the sequence name
		if (seqNames[i].length() > 0) {
			x_coord = x - HEATMAP_GAP_BW_TEXT;
			y_coord = y + block_size*0.8;
			outputSeqName(fout, x_coord, y_coord, 0, seqNames[i], HEATMAP_FONT, font_size, HEATMAP_TEXT_ALIGN);
		}

		for (j=0; j<i; j++) {
			// print out the squares representing the Cij values
			outputSquare(fout, x, y, getColor(cij[i*seqNum+j]), block_size);
			x += block_size;
		}

		// print out the gray square
		outputSquare(fout, x, y, HEATMAP_GRAY_COLOR, block_size);
		x += block_size;
		
		for (j=i+1; j<seqNum; j++) {
			// print out the squares representing the Cij values
			outputSquare(fout, x, y, getColor(cij[i*seqNum+j]), block_size);
			x += block_size;
		}
		
		y += block_size;
	}
	
	// output the legend
	y = HEATMAP_TOP_BORDER;
	x = HEATMAP_LEFT_BORDER + longestTextLen*font_size/2.0 + HEATMAP_GAP_BW_TEXT + 
				block_size*seqNum + HEATMAP_GAP_BW_LEGEND;
	outputLegend(fout, x, y);
	
	// output the tail
	outputSVGTail(fout);

	fout.close();
}
