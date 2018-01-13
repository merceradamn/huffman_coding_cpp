/*
Name:		Adam Mercer
Project:	Huffman Algorithm Implementation
Descript:	Using a body of text, implement the Huffman's Algorithm to create 
			a compressed version of the text.

			Also added is a way of not needing a separate file to decompress the
			file since the prefix table is added to the file along with the compressed
			text. The only hang up is that it becomes somewhat ridiculous to try and
			compress a file with only several bytes of information.

			The first step gets the file and counts the occurence of every character
			and stores that in a map. We then use that map and flip it so the frequency
			maps to the character.

Resources:	http://www.geeksforgeeks.org/greedy-algorithms-set-3-huffman-coding/
			This code is contributed by Aditya Goel and adapted by Adam Mercer
			Structure for building a Huffman Tree - Adapted to handle using a vector

Notes:		
*/

#include <iostream>				// For using cout
#include <string>				// String-handling header
#include <iomanip>				// Formatting header
#include <cmath>				// Math functions header
#include <fstream>				// File-handling header
#include <queue>				// Queue container header
#include <vector>				// Vector container header
#include <map>					// Map container header
#include <stack>				// Stack container header
#include <bitset>
using namespace std;

// A Huffman tree node
struct MinHeapNode {
	char data;											// One of the input characters
	unsigned freq;										// Frequency of the character
	MinHeapNode *left, *right;							// Left and right child

	MinHeapNode(char data, unsigned freq) {
		left = right = NULL;
		this->data = data;
		this->freq = freq;
	}
};

// For comparison of two heap nodes (needed in min heap)
struct compare {
	bool operator()(MinHeapNode* l, MinHeapNode* r) {
		return (l->freq > r->freq);
	}
};

map<char, int> fileReader(map<char, int> &map, const string &filename);
void printCodes(struct MinHeapNode* root, string str, map<char, string> &sM);
priority_queue<MinHeapNode*, vector<MinHeapNode*>, compare> HuffmanCodesv2(vector<pair<char, int>> chFq);
void encode(map<char, string> stMp, string &filename, string &textenc);
void decodeStandalone(const string &textenc);

int main(){
/*************************************************************************************************/
/*	VARIABLES AND DECLARATIONS																	 */
/*************************************************************************************************/				
	// Default filenames for testing purposes
	string filenames[7];
	filenames[0] = "huffmantext.txt";						// Name of file 1
	filenames[1] = "Alice’s Adventures in Wonderland.txt";	// Name of file 2
	filenames[2] = "Pride and Prejudice.txt";				// Name of file 3
	filenames[3] = "A Tale of Two Cities.txt";				// Name of file 4
	filenames[4] = "The Adventures of Tom Sawyer.txt";		// Name of file 5
	filenames[5] = "huffmanAlg_cpp.txt";					// Name of file 6
	//filenames[6] = "Testing.txt";

	int x = 0;												// Choice of file
	bool good = false;										// Bool for good or bad file choice

	// Related to the encode/decode process
	string encoded, decoded;								// Filename variables for encode/decode
	map<char, int> fileFreq;								// Map holding our frequency chart
	map<char, int>::iterator fiFr;							// Iterator for fileFreq map
	vector<pair<char, int>> charFreq;						// Vector of char/freq pairs
	priority_queue<MinHeapNode*, vector<MinHeapNode*>, compare> minHeap;
	map<char, string> strMap;								// Char/Var-length code map
	bool redo = false;

/*************************************************************************************************/
/*	HUFFMAN TREE BUILDING																		 */
/*************************************************************************************************/
	int f = 9;	
	while (!good) {
		// Check which file we want to encode
		if (f == 9) {
			cout << "Files to test encode/decode:\n";
			for (int i = 1; i < 8; i++)
				cout << setw(36) << filenames[i - 1] << " (" << i << ")\n";
			f = 0;
		}
		cout << "What file do you want to encode/decode(1-7): ";
		cin >> x;
		if (x > 0 && x <= 7) {
			good = true;
			x--;		// Set the choice of file to encode/decode
		}
		else if (x == 0) {
			cout << "Ending encode/decode program.\n";
			system("pause");
			return 0;
		}
		else if (x == 9) {
			f = x;		// Reprint the file names
			continue;
		}
		else if (x < 0 || x > 7) {
			cout << "Error. Choice of file.\n";
			good = false;
			continue;
		}

		// Read the file and get the frequency of each character
		fileReader(fileFreq, filenames[x]);

		// Check if the file is empty
		if (fileFreq.empty()) {
			cout << "File is empty. Can't encode an empty file.\n";
			good = false;
			continue;
		}

		// Add the map to a vector to preserve integrity
		// Having a map with characters that map to similar values
		// seems dangerous enough to try to prevent this with a vector
		fiFr = fileFreq.begin();
		for (; fiFr != fileFreq.end(); fiFr++) {
			charFreq.push_back(make_pair(fiFr->first, fiFr->second));
		}

		// Passes the vector we made previously and builds the tree using
		// the struct object we declared further above and then store that
		// in a variable we declared in main
		minHeap = HuffmanCodesv2(charFreq);					// Build the heap holding the char/freqs

		// Check the status of the heap, if we've got a size to move on
		if (minHeap.empty()) {
			cout << "File can't be encoded. Returning to prompt.\n";
			fileFreq.clear();									// Clear the frequency chart
			charFreq.clear();									// Clear vector holding char/freq pairs
			good = false;
			continue;
		}

		// Check if the heap
		printCodes(minHeap.top(), "", strMap);				// Store the heap data in a map for quick search

		// Set the name of what we want encoded to be called
		encoded.assign(filenames[x], 0, filenames[x].size() - 4);
		encoded += " - encoded.bin";
		cout << "Name of encoded file is " << encoded << endl;

		encode(strMap, filenames[x], encoded);				// Tree, Orig File, Encoded Filename
		decodeStandalone(encoded);							// Pass name of encoded file

		// Clear the variables so we can encode and decode again
		fileFreq.clear();									// Clear the frequency chart
		charFreq.clear();									// Clear vector holding char/freq pairs
		strMap.clear();										// Map holds char & bit code
		minHeap.pop();										// Clear the heap holding our char/freqs
		good = false;										// Clear the check for file

	}
	return 0;
}

/*************************************************************************************************/
/*	FUNCTIONS																					 */
/*************************************************************************************************/

// Reads a file and counts the frequency of each unique character
// Stores them in a map
map<char, int> fileReader(map<char, int> &fF, const string &filename){

	// Try to set the locale to UTF8, if it fails just use whatever is default
	locale loc;
	try {
		loc = std::locale("en_US.UTF8");
	}
	catch (std::runtime_error) {
		loc = std::locale(loc, "", std::locale::ctype);
	}

	map<char, int>::iterator it = fF.begin();			// Iterator to the beginning of map
	char ch;											// Holds a character for processing
	ifstream fileIn;									// Create an ifstream object
	fileIn.open(filename);								// Open the file

	if (!fileIn.is_open()){								// Make sure the file is opened
		cout << "File not opened." << endl;
		return fF;										// Return out of the function if not
	}

	// Start getting a character at a time from the file
	while (fileIn.get(ch)){								// Get the next character in the file
		it = fF.find(ch);								// Check the map for the character
		if (it != fF.end()) {							// If found the character
			fF[ch]++;									// Increment it's frequency
		}
		else {
			fF[ch] = 1;									// Add character to map and set freq to 1
		}
	}

	fileIn.close();
	return fF;
}

// Prints huffman codes from the root of Huffman Tree
// Also stores the data in a map for easy access
void printCodes(struct MinHeapNode* root, string str, map<char, string> &sM) {
	// Should return if we don't have a tree
	if (!root)
		return;

	// $ is used in a parent node to indicate that we don't have a valid code yet
	if (root->data != '$') {
		//cout << root->data << ": " << setw(10) << str << "\t";
		//cout << root->freq << "\n";
		// Store the data in a map
		sM.insert(pair<char, string>(root->data, str));
	}

	// Recursively call this function to print each code for it's character
	// Passes the next node down and either "0" or "1" based on direction taken
	printCodes(root->left, str + "0", sM);
	printCodes(root->right, str + "1", sM);
}

// This variation of the Huffman Tree building uses a vector of pairs to build the tree
// Rather than print the code afterwards and not do anything with the queue after printing
// this version will return the queue and skip the print step for later
priority_queue<MinHeapNode*, vector<MinHeapNode*>, compare> HuffmanCodesv2(vector<pair<char, int>> chFq) {
	// Variables for building the tree
	vector<pair<char, int>>::iterator chFqIt = chFq.begin();
	struct MinHeapNode *left, *right, *top;

	// Create a min heap & inserts all characters of chFq
	// chFq is the vector of pairs with character/int as the members
	priority_queue<MinHeapNode*, vector<MinHeapNode*>, compare> minHeap;

	if (chFq.size() > 1) {
		for (unsigned int i = 0; i < chFq.size(); ++i) {
			minHeap.push(new MinHeapNode(chFqIt->first, chFqIt->second));
			chFqIt++;
		}

	}
	else if (chFq.size() == 1) {
		cout << "File too small to get any compression. Returning to prompt.\n";
		return minHeap;
	}

	// Iterate while size of heap isn't equal to 1
	while (minHeap.size() != 1) {
		// Extract the two minimum freq items from min heap
		left = minHeap.top();
		minHeap.pop();

		right = minHeap.top();
		minHeap.pop();

		// Create a new internal node with frequency equal to the
		// sum of the two nodes frequencies. Make the two extracted
		// node as left and right children of this new node. Add
		// this node to the min heap
		// '$' is a special value for internal nodes, not used
		top = new MinHeapNode('$', left->freq + right->freq);
		top->left = left;
		top->right = right;
		minHeap.push(top);
	}

	// Print Huffman codes using the Huffman tree built above
	//printCodes(minHeap.top(), "");

	return minHeap;
}

// Encode the original text and barring any failed character reads we should be able to get
// out an encoding of the original file in a much leaner version than before
void encode(map<char, string> stMp, string &filename, string &textenc) {

	// Open the encode file for output
	fstream encode;										// Create our fstream object
	encode.open(textenc, ios::out | ios::binary | ios::trunc);	// Open the file and trunc any data
	queue<char> encodedBytes;							// Holds a queue of our bytes to write

	// Open the original text
	ifstream fileIn;									// Create an ifstream object
	fileIn.open(filename);								// Open the file

	if (!fileIn.is_open()) {							// Make sure the file is opened
		cout << "File not opened." << endl;
		cout << "Encoding Halted!" << endl;
		return;											// Return out of the function if not
	}

	// Encode each character from the input file
	char ch; char temp = 0;								// Character variables
	unsigned int numToCh = 0;							// Running number for byte crunching
	unsigned int charCount = 0;							// Holds the number of characters encoded
	unsigned int bitCount = 0;
	queue<char> bits;									// Queue for holding our bits
	map<char, string>::iterator stMpIT;					// Iterator for the string map
	while (fileIn.get(ch)) {							// Get the next character in the file
		stMpIT = stMp.find(ch);							// Use iterator to try to find character
		if (stMpIT != stMp.end()) {						// If iterator not at end, process char
			charCount++;								// Increment our char count for header info
			
			// Store the bits into a queue for byte crunching
			bitCount += stMpIT->second.size();			// Add the number of bits we got

			for (auto & ch : stMpIT->second) {			// For each character in the string
				bits.push(ch);							// Push a bit on the queue
			}
			
			// Check the queue for enough bits and handle processing a byte
			if (bits.size() >= 8) {						// Check if we have 8 bits
				for (int i = 7; i > -1; i--) {			// If 8 bits then store in char and write
					// Pop off 8 bits
					temp = bits.front();				// Store the bit in a variable
					bits.pop();							// Pop off the front bit
					if (temp == '1') {					// Calculate if bit is 1
						numToCh += pow(2, i);			// Add to the running total in nToC
					}
				}

				// Write the byte to the file
				encodedBytes.push(numToCh);				// Queue the byte
				numToCh = 0;							// Reset the numtochar variable
			}
		}
		else {
			cout << "Invalid character found! => ";
			cout << "{" << ch << "}" << endl;
		}
	}

	// Finished reading the file so we need to check for leftover bits
	// Also we need to check to see if we never had enough bits to fill a byte
	if (bits.size() > 0) {
		int buffbits = 7 - bits.size();						// Check if we still have bits queue'd
		for (int i = buffbits; i > -1; i--) {				// Format the remaining bits with buffer
			if (!bits.empty()) {
				temp = bits.front();							// Store the bit in a variable
				bits.pop();										// Pop off the front bit
				if (temp == '1') {								// Calculate if bit is 1
					numToCh += pow(2, i);						// Add to the running total in nToC
				}
			}
		}

		// Write the last byte to the file after converting
		encodedBytes.push(numToCh);							// Queue the byte
	}

	// Use an iterator to write the tree information to the file for decoding
	encode << stMp.size() << "\t" << bitCount << endl;			/*EDIT: encodedBytes.size()*/
	for (stMpIT = stMp.begin(); stMpIT != stMp.end(); stMpIT++) {
		encode << stMpIT->first << stMpIT->second << " ";
	}
	encode << endl;

	// Write the contents of the encoded byte queue now
	while (!encodedBytes.empty()) {						// Check for an empty enc byte queue
		encode << encodedBytes.front();					// Write the char to the file
		encodedBytes.pop();								// Pop the character off the queue
	}

	// Close the files
	fileIn.close();
	encode.close();
}

// Decode an encoded text file
// Uses nothing but the file as we wrote the code tree to the file
// beforehand so we just need to process the file
void decodeStandalone(const string &textenc) {

	// Variable list
	queue<char> bitQueue;								// Holds all the bits after decoding
	map<string, char> codeMap;							// Holds the map for char/code pairs
	map<string, char>::iterator CMIT;					// codeMap iterator
	char readCh;										// Character read from encoded file
	string chCode;										// Bit code for the character
	int codeCount;										// Count of how many Huffman Codes to read
	int bitCount;										// Number of bytes to read and decode
	stack<int> bits;									// Stack holding bits for reversing							
	char nextCh;										// Holds a character to read in

	// Open the file we encoded
	fstream readEncode;									// Create fstream object
	readEncode.open(textenc, ios::in | ios::binary);	// Open the file in read mode

	if (!readEncode.is_open()) {						// Make sure the file is opened
		cout << "File not opened." << endl;
		cout << "Decoding Halted!" << endl;
		return;											// Return out of the function if not
	}

	// Open the file that we're going to write the decoded text to
	fstream writeDecode;								// Decoded text object
	string decodedText;									// Filename for decoded text
	decodedText.assign(textenc, 0, textenc.size() - 14);
	decodedText += " - decoded.txt";
	cout << "Name of decoded file is " << decodedText << endl;
	writeDecode.open(decodedText, fstream::out, fstream::trunc);

	// Check if file has data, if so get code/byte counts and prepare to get code table
	if (!readEncode.eof()) {							// Make sure the file isn't empty
		readEncode >> codeCount;						// Get the number of Huffman codes for the file
		readEncode >> bitCount;						// Get the number of bytes we need to convert
		readEncode.get(readCh);							// Get the new line character
	}
	else {
		cout << "Empty File. Decode failed." << endl;	// Scrub decode if we don't have valid file
		return;
	}

	// Get the Huffman code for each character so we can decode the file
	for (int i = 0; i < codeCount; i++) {
		readEncode.get(readCh);							// Get the character
		readEncode >> chCode;							// Get the bit code
		readEncode.get();								// Get the trash space character
		codeMap.insert(pair<string, char>(chCode, readCh));
	}

	// Start getting characters from the file and converting them to binary
	int chToDec = 0;									// Integer for converting char to binary
	int counter = bitCount;
	readEncode.get();									// Discard the newline after char/code table
	while (readEncode.get(nextCh)) {					// Read a character if it exists
		int numBits;
		bitset<8> chToBits(nextCh);						// Get the related bits of a character

		if (counter >= 8)								// If the bit counter is bigger than 8	
			numBits = 8;								// Set the bits so we pull all 8
		else if (counter <= 4)							// If the bit counter is smaller than 8
			numBits = 4;								// Set the bits to 4, only get 4 bits

		if (numBits == 8) {								// If we need all 8 bits
			for (int i = 7; i > -1; i--) {				// Store the bits to queue, in the right order
				bitQueue.push(chToBits[i]);
				counter--;
			}
		}
		else if (numBits == 4) {						// If we just need 4 bits
			for (int i = 3; i > -1; i--) {				// Store the bits to queue, in the right order
				bitQueue.push(chToBits[i]);
				counter--;
			}
		}
	}

	// Get the first 3 bits from the queue and check the map
	// If no match add another bit
	// When you get a match pop off the usable bits from the queue
	string bitsToCheck;									// String to compare to map values
	CMIT = codeMap.begin();								// Set the iterator for comparing to the map
	while (!bitQueue.empty() || bitCount > -1) {
		CMIT = codeMap.find(bitsToCheck);				// Look for a matching key
		if (CMIT == codeMap.end()) {					// If no key found
			if (bitQueue.empty() && !bitsToCheck.empty()){// Check for more bits on the queue
				cout << "Ran out of bits to read." << endl;
				return;
			}
			else if(!bitQueue.empty()){					// Check if we have bits to queue
				// Gotta adjust the bits to "0"/"1" in ascii for comparing
				bitsToCheck += bitQueue.front() + 48;	// Add the character to the string
				bitQueue.pop();							// Pop the bit we just got off the queue
				bitCount--;
			}
			else if (bitQueue.empty() && bitsToCheck == "") {// Check if queue and check are empty
				writeDecode.close();					// Close the file
				return;									// Return since we finished
			}
		}
		else if (CMIT->first == bitsToCheck) {			// We got a matching character from the map
			bitsToCheck.clear();						// Clear the string
			writeDecode << CMIT->second;				// Write the character to decoded file
		}
	}

}