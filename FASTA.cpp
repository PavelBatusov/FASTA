#include <map>
#include <time.h> 
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "sqlite3.h"

#define DB_CONFIG_FILE_NAME "DB_config.txt"
#define DB_NAME "FASTA.db"

bool filter1(int d_count, int m_count, 
						 const std::string &seq1, const std::string &seq2) 
{
	//alloc=======================================================================
	int n = seq1.length(), m = seq2.length();
	int* matrix = new int[n * m * sizeof(int)];
	int maxlen = 0, count = 0;
	//initialize==================================================================
	for (int i = 0; i < seq1.length(); i++) matrix[i*m+0] = (seq1[i] == seq2[0]);
	for (int j = 0; j < seq2.length(); j++) matrix[0*m+j] = (seq1[0] == seq2[j]);
	//calc========================================================================
	for (int i = 1; i < seq1.length(); i++) {
		for (int j = 1; j < seq2.length(); j++) {
			if (seq1[i] != seq2[j]) {
				if (matrix[(i-1)*m+j-1] > 1) count++;
			} else {
				matrix[i*m+j] = matrix[(i-1)*m+j-1] + 1;
				if (matrix[i*m+j] > maxlen) maxlen = matrix[i*m+j];
			}
		}
	}
	//free========================================================================
	delete[] matrix;
	//check=======================================================================
	return (count < d_count || maxlen < m_count);
}

static int callback(void* NotUsed, int argc, char** argv, char** azColName){
	for (int i = 0; i < argc; i++) {
		printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
	}
	printf("\n");
	return 0;
}

void InsertDB(char* file_name, int magic_num) {
	//open DB=====================================================================
	sqlite3 *db;
	int rc = sqlite3_open(DB_NAME, &db);
	if( rc ) {
		printf("Can't open database: %s\n", sqlite3_errmsg(db));
		sqlite3_close(db);
		return;
	}
	//open FileName===============================================================
	std::ifstream input_seq(file_name);
	std::string name = "";
	std::string seq = "";
	char c;
	//insert input sequences------------------------------------------------------
	while (input_seq.get(c)) {
		switch (c) {
			case '>':
				if (seq.length() && name.length()) {
					//insert seq---------------------------------------------------------- 
					sqlite3_stmt *stmt;
					const char *pzTest;
					std::string insert("INSERT INTO Sequence (description, string) VALUES\
(?, ?)");
					rc = sqlite3_prepare_v2(db, insert.c_str(), insert.length(), 
																	&stmt, &pzTest);
					if( rc == SQLITE_OK ){
						sqlite3_bind_text(stmt, 1, name.c_str(), name.length(), 0);
						sqlite3_bind_text(stmt, 2, seq.c_str(), seq.length(), 0);
						sqlite3_step(stmt);
						sqlite3_finalize(stmt);
					} else { 
						printf("Can't prepare sqlite3_prepare_v2 'INSERT INTO Sequence (des\
cription, string) VALUES(%s, %s)'\n", name.c_str(), seq.c_str()); 
					}
					//inser subseq--------------------------------------------------------
					int mainID = sqlite3_last_insert_rowid(db);
					for (int i = 0; i < seq.length() - magic_num + 1; i++) {
						insert = "INSERT INTO SubSequence (string, position, baseString) VA\
LUES(?, ?, ?)";
						rc = sqlite3_prepare_v2(db, insert.c_str(), insert.length(), 
																		&stmt, &pzTest);
						if( rc == SQLITE_OK ){
							sqlite3_bind_text(stmt, 1, seq.substr(i, magic_num).c_str(), 
																magic_num, 0);
							sqlite3_bind_int(stmt, 2, i);
							sqlite3_bind_int(stmt, 3, mainID);
							sqlite3_step(stmt);
							sqlite3_finalize(stmt);
						} else { 
							printf("Can't prepare sqlite3_prepare_v2 'INSERT INTO SubSequence\
 (string, position, baseString) VALUES(%s, %d, %d)'\n", 
										 seq.substr(i, magic_num).c_str(), i, mainID); 
						}
					}
				}
				std::getline(input_seq, name);
				seq = "";
				break;
			case '\t': 
			case '\n': 
			case '\r': 
			case  ' ': 
				//skip ws---------------------------------------------------------------
				break;
			default:
				seq += c;
		}
	}
	//close file==================================================================
	input_seq.close();
	//close DB====================================================================
	sqlite3_close(db);
}

void CreateDB(char* file_name, char* magic_num) {
	//variables to work with database=============================================
	sqlite3 *db;
	char *zErrMsg = 0;
	int rc;
	std::string dbName(DB_NAME);
	std::ifstream fin("FASTAquery.sql");
	//loading SQL query to create DB----------------------------------------------
	std::string query;
	std::getline(fin, query, '$');
	//open && create DB-----------------------------------------------------------
	rc = sqlite3_open(dbName.c_str(), &db);
	if( rc ) {
		printf("Can't open database: %s\n", sqlite3_errmsg(db));
		sqlite3_close(db);
		return;
	}
	rc = sqlite3_exec(db, query.c_str(), callback, 0, &zErrMsg);
	if( rc != SQLITE_OK ) {
		printf("SQL error: %s\n", zErrMsg);
		sqlite3_free(zErrMsg);
		sqlite3_close(db);
		return;
	}
	sqlite3_close(db);
	//insert data to base=========================================================
	InsertDB(file_name, atoi(magic_num));
	//creating DB config file=====================================================
	std::ofstream conf_file(DB_CONFIG_FILE_NAME);
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	conf_file << magic_num << '\n';
	conf_file << "Creation time: " << asctime(timeinfo);
	conf_file.close();
}

void AddDB(char* file_name) {
	//trying to open DB config file
	std::ifstream conf_file(DB_CONFIG_FILE_NAME);
	int magic_num = 0;
	conf_file >> magic_num;
	conf_file.close();
	//adding sequences to DB
	InsertDB(file_name, magic_num);
}

void Search(char* seq_file_name) {
	//read input seq==============================================================
	std::ifstream input_seq(seq_file_name);
	std::string seq = "";
	//read seq--------------------------------------------------------------------
	std::getline(input_seq, seq);
	//read extra params-----------------------------------------------------------
	int d_count, m_count;
	input_seq >> d_count >> m_count;
	//read alphabet && score matrix-----------------------------------------------
	//....
	//close file------------------------------------------------------------------
	input_seq.close();
	//open DB config file=========================================================
	std::ifstream conf_file(DB_CONFIG_FILE_NAME);
	int magic_num = 0;
	conf_file >> magic_num;
	conf_file.close();
	//start work with DB==========================================================
	sqlite3 *db;
	int rc = sqlite3_open(DB_NAME, &db);
	if( rc ){
		printf("Can't open database: %s\n", sqlite3_errmsg(db));
		sqlite3_close(db);
		return;
	}
	//for substr------------------------------------------------------------------
	//hashmap for all strings from base who is good
	std::map<int, std::string> dump;
	for (int i = 0; i < seq.length() - magic_num + 1; i++) {
		std::string substr = seq.substr(i, magic_num);
		std::string select = "SELECT ID, string FROM sequence INNER JOIN (SELECT DI\
STINCT baseString FROM SubSequence WHERE string = ? ) AS keys  ON sequence.id =\
 keys.baseString";
		sqlite3_stmt *stmt;
		const char *pzTest;
		rc = sqlite3_prepare_v2(db, select.c_str(), select.length(), &stmt,&pzTest);
		if( rc == SQLITE_OK ) {
			sqlite3_bind_text(stmt, 1, substr.c_str(), substr.length(), 0);
			while (sqlite3_step(stmt) != SQLITE_DONE) {
				int id_seq2 = sqlite3_column_int(stmt, 0); 
				std::string seq2((char*)sqlite3_column_text(stmt, 1));
				dump.insert(std::pair<int, std::string>(id_seq2, seq2.c_str()));
			}
			sqlite3_finalize(stmt);
		} else {
			//some shit happened
			printf("Can't prepare sqlite3_prepare_v2 'SELECT ID, string FROM sequence\
 INNER JOIN (SELECT DISTINCT baseString FROM SubSequence WHERE string = %s) AS \
keys  ON sequence. id = keys.baseString'", substr.c_str());
		}
	}
	//start filtering============================================================= 
	printf("loading %lu sequences from database\n", dump.size());
	//first filter - diagonals count && max length--------------------------------
	for (auto it = dump.begin(); it != dump.end(); ) {
		if (filter1(d_count, m_count, seq, (*it).second)) {
			it = dump.erase(it);
		} else it++;
	}
	printf("after first filter: %lu sequences\n", dump.size());
	//second filter---------------------------------------------------------------
	//...
	//S-W for winners
	//...
}

int main(int argc, char** argv) {
	if (argc < 2) {
		printf("Too few arguments\n");
		printf("\tpossible arguments:\n");
		printf("\t\tc <file name> <magic num> - create new DB from \
<file name>with <magic num>\n");
		printf("\t\ta <file name> - add sequences from <file name> to \
exsistingDB\n");
		printf("\t\ts <file name> - search seq from <file name> in \
exsisting DB\n");
		return 0;
	}
	switch (*argv[1]) {
		case 'c':
			//create
			if (argc > 3) CreateDB(argv[2], argv[3]);
			else printf("Not enough arguments\n");
			break;
		case 'a':
			//add
			if (argc > 2) AddDB(argv[2]);
			else printf("Not enough arguments\n");
			break;
		case 's':
			//search
			if (argc > 2) Search(argv[2]);
			else printf("Not enough arguments\n");
			break;
		default:
			printf("Unknown command %s\n", argv[1]);
	}
	return 0;
}