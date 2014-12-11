#include <map>
#include <time.h> 
#include <string>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "sqlite3.h"

#define DB_CONFIG_FILE_NAME "DB_config.txt"
#define DB_NAME "FASTA.db"

//результаты прохода Смита-Ватермана
struct SWres {
	int score;
	int i, j;
	int* way;
	std::string seq1, seq2, name;
	
	SWres(int scr, int pos_i, int pos_j, int* w, std::string s1, std::string s2):
		score(scr), i(pos_i), j(pos_j), way(w), seq1(s1), seq2(s2) {};	
	
	~SWres() {}
	
	bool operator < (const SWres &other) const { return (score > other.score); }
	
	void FreeSharedMem() { if (way) delete[] way; }
};

//удаление пробельных символов
void SpaceErase(std::string &str) {
	str.erase(std::remove_if(str.begin(), str.end(),
	[](char c){ return (c == ' ' || c == '\n' || c == '\t'); }), str.end());
} 

//чтение алфавита, матрицы очков и штрафа за геп 
void GetScoreMatrix(std::string &alphabet, int* index_arr, int* &score_matrix, 
										int &penalty, std::istream& fs) 
{
	fs >> std::ws;
	
	std::getline(fs, alphabet);
	SpaceErase(alphabet);
	
	//факторизация алфавита
	for (int i = 0; i < alphabet.length(); i++) {
		index_arr[alphabet[i]] = i;
	}
	
	score_matrix = new int [alphabet.length() * alphabet.length()];
	for (int i = 0; i < alphabet.length(); i++) {
		for (int j = 0; j < alphabet.length(); j++) {
			fs >> score_matrix[i * alphabet.length() + j];
		}
	}
	
	fs >> penalty;
}

bool filter1(int d_count, int m_count, 
						 const std::string &seq1, const std::string &seq2) 
{
	//alloc=======================================================================
	int n = seq1.length(), m = seq2.length();
	int* matrix = new int[n * m * sizeof(int)];
	int maxlen = 0, count = 0;
	
	//инициализация граничных строк===============================================
	for (int i = 0; i < seq1.length(); i++) 
		matrix[i*m+0] = (seq1[i] == seq2[0]);
	
	for (int j = 0; j < seq2.length(); j++) 
		matrix[0*m+j] = (seq1[0] == seq2[j]);
	
	//calc========================================================================
	for (int i = 1; i < seq1.length(); i++) {
		for (int j = 1; j < seq2.length(); j++) {
			if (seq1[i] != seq2[j]) {
				if (matrix[(i-1)*m + j-1] > 1) count++;
				matrix[i*m + j] = 0;
			} else {
				matrix[i*m + j] = matrix[(i-1)*m + j-1] + 1;
				if (matrix[i*m + j] > maxlen) maxlen = matrix[i*m + j];
			}
		}
	}
	
	//free========================================================================
	delete[] matrix;
	
	//check=======================================================================
	return (count < d_count || maxlen < m_count);
}

bool filter2(int score_infinum, int* score_matrix,
						 int* index_arr, int alphabetLength, 
						 const std::string &seq1, const std::string &seq2) 
{
	//alloc=======================================================================
	int n = seq1.length(), m = seq2.length();
	int* matrix = new int[n * m * sizeof(int)];
	int maxlen = 0, count = 0;
	
	//инициализация границ========================================================
	for (int i = 0; i < seq1.length(); i++) {
		int scr = score_matrix[index_arr[seq1[i]] * alphabetLength +
													 index_arr[seq1[0]]
													];
		matrix[i*m+0] = (scr > 0) ? scr : 0;
	}
	
	for (int i = 0; i < seq2.length(); i++){
		int scr = score_matrix[index_arr[seq1[0]] * alphabetLength + 
																 index_arr[seq1[i]]
																];	
		matrix[0*m+i] = (scr > 0) ? scr : 0;
	}
	
	//calc========================================================================
	int maxScore = 0;
	for (int i = 1; i < seq1.length(); i++) {
		for (int j = 1; j < seq2.length(); j++) {
			if (seq1[i] == seq2[j]) {
				matrix[i*m + j] =
					score_matrix[index_arr[seq1[i]] * alphabetLength + 
											 index_arr[seq2[j]]
											] +
					matrix[(i-1)*m + j-1];
				if (matrix[i*m+j] > maxScore)
					maxScore = matrix[i*m+j];
			} else 
				matrix[i*m + j] = 0;
		}
	}
	
	//free========================================================================
	delete[] matrix;
	
	//check=======================================================================
	return (maxScore < score_infinum);
}

//колбек для работы с базой
static int callback(void* NotUsed, int argc, char** argv, char** azColName){
	for (int i = 0; i < argc; i++) {
		printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
	}
	printf("\n");
	return 0;
}

void GetAllign(int* way, int i, int j, 
							 const std::string &seq1, const std::string &seq2, 
							 std::string &result1, std::string &result2) 
{
	int m = seq2.length();
	
	//добавление правых "хвостов"
	if (j == m) result1 = seq1.substr(i, seq1.length() - i);
	else result2 = seq2.substr(j, seq2.length() - j);
	
	//обратный ход по таблице
	while (i && j) {
		if (way[i*m + j] == 1) { 
			result1.insert(0, 1, seq1[--i]); 
			result2.insert(0, 1, seq2[--j]);  
		} else if (way[i*m + j] == 2) { 
			result1.insert(0, 1, '-'); 
			result2.insert(0, 1, seq2[--j]);  
		} else { 
			result1.insert(0, 1, seq1[--i]); 
			result2.insert(0, 1, '-');  
		}
	}
	
	while (i) {
		result1.insert(0, 1, seq1[--i]);
		result2.insert(0, 1, ' ');
	}
	
	while (j) {
		result1.insert(0, 1, ' ');
		result2.insert(0, 1, seq2[--j]);		
	}
}

SWres SmithWaterman(int* score_matrix, int* index_array, int alpha_len,
										int penalty, const std::string &seq1, const std::string &seq2) 
{
	int n = seq1.length(), m = seq2.length();
	int* score = new int [(n + 1) * (m + 1)];
	int* way = new int [(n + 1) * (m + 1)];
	
	for (int i = 0; i <= n; i++) {
		score[i * m] = 0;
	}
	
	for (int j = 0; j <= m; j++) {
		score[j] = 0;
	}
	
	for (int i = 1; i < n + 1; i++) {
		for (int j = 1; j < m + 1; j++) {
			//оставить все как есть
			int way1 = score[(i-1)*m + j-1] + 
				score_matrix[index_array[seq1[i-1]]*alpha_len + index_array[seq2[j-1]]];
			//порвать seq1
			int way2 = score[i*m + j-1] + penalty;
			//порвать seq2
			int way3 = score[(i-1)*m + j] + penalty;
			//выбираем максимум
			if (way1 >= way2 && way1 >= way3) {
				score[i*m + j] = way1;
				way[i*m + j] = 1;
			} else if (way2 >= way1 && way2 >= way3) {
				score[i*m + j] = way2;
				way[i*m + j] = 2;
			} else {
				score[i*m + j] = way3;
				way[i*m + j] = 3;
			}
		}
	}
	
	int max_score = score[m], fi = 0, fj = 0;
	
	for (int i = 1; i <= n; i++) {
		if (score[m*i + m] > max_score) {
			max_score = score[m*i + m];
			fi = i;
		}
	}
	
	for (int j = 1; j <= m; j++) {
		if (score[m*n + j] > max_score) {
			max_score = score[m*n + j];
			fj = j;
		}
	}
	
	return (score[m*fi + m] > score[m*n + fj]) ? 
		SWres(max_score, fi, m, way, seq1, seq2) :
		SWres(max_score, n, fj, way, seq1, seq2);
}

void InsertDB(char* file_name, int magic_num) {
	//открываем базу==============================================================
	sqlite3 *db;
	int rc = sqlite3_open(DB_NAME, &db);
	if( rc ) {
		printf("Can't open database: %s\n", sqlite3_errmsg(db));
		sqlite3_close(db);
		return;
	}
	
	//открываем файл с последовательностями=======================================
	std::ifstream input_seq(file_name);
	std::string name = "";
	std::string seq = "";
	char c;
	
	//вставка последовательностей в базу------------------------------------------
	while (input_seq.get(c)) {
		switch (c) {
			case '>':
				//начинается новая последовательность-----------------------------------
				if (seq.length() && name.length()) { //если уже что-то считано
					//вставить последовательность в базу
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
					
					//вставка подпоследовательностей в базу
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
			//опускаем пробельные символы---------------------------------------------
			case '\t': 
			case '\n': 
			case '\r': 
			case  ' ': 
				break;
			//чтение последовательности-----------------------------------------------
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
	//переменные для работы с базой===============================================
	sqlite3 *db;
	char *zErrMsg = 0;
	int rc;
	std::string dbName(DB_NAME);
	std::ifstream fin("FASTAquery.sql");
	
	//загружаем SQL запрос на создание базы---------------------------------------
	std::string query;
	std::getline(fin, query, '$');
	
	//открываем и создаем базу----------------------------------------------------
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
	
	//заносим данные в базу=======================================================
	InsertDB(file_name, atoi(magic_num));
	
	//создание конфиг. файла базы=================================================
	std::ofstream conf_file(DB_CONFIG_FILE_NAME);
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	conf_file << magic_num << '\n';
	conf_file << "Creation time: " << asctime(timeinfo) << std::endl;
	conf_file.close();
}

void AddDB(char* file_name) {
	//пытаемся открыть конфиг. файл базы
	std::ifstream conf_file(DB_CONFIG_FILE_NAME);
	int magic_num = 0;
	conf_file >> magic_num;
	conf_file.close();
	
	//вставляем данные в базу
	InsertDB(file_name, magic_num);
	
	//немного логов :)
	std::fstream fs(DB_CONFIG_FILE_NAME, std::ios::out|std::ios::app);
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	fs << "Modifed: " << asctime(timeinfo) << std::endl;
	fs.close();
}

void Search(char* seq_file_name) {
	//входная последовательность==================================================
	std::ifstream input_seq(seq_file_name);
	std::string seq = "";
	std::getline(input_seq, seq);
	
	//параметры поиска------------------------------------------------------------
	int d_count, m_count, p_count;
	input_seq >> d_count >> m_count >> p_count;
	
	//чтение алфавита и матрицы очков---------------------------------------------
	std::string alphabet;
	int index_arr[128];
	int* score_matrix = NULL;
	int penalty;
	
	GetScoreMatrix(alphabet, index_arr, score_matrix, penalty, input_seq);
	
	//чтение закончено------------------------------------------------------------
	input_seq.close();
	
	//инициализация работы с базой================================================
	//считываем информацию из конфиг. файла
	std::ifstream conf_file(DB_CONFIG_FILE_NAME);
	int magic_num = 0;
	conf_file >> magic_num;
	conf_file.close();
	
	//открытие базы===============================================================
	sqlite3 *db;
	int rc = sqlite3_open(DB_NAME, &db);
	if( rc ){
		printf("Can't open database: %s\n", sqlite3_errmsg(db));
		sqlite3_close(db);
		return;
	}
	
	//для всех подстрок входной последовательности--------------------------------
	//составляем словарь со строчками из базы, имеющии такие подпоследовательности
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
	
	//фильтрация первичной выборки================================================ 
	printf("loading %lu sequences from database\n", dump.size());
	
	//запуск первого фильтра - проверка на число диагоналей и макс. длинну--------
	for (auto it = dump.begin(); it != dump.end(); ) {
		if (filter1(d_count, m_count, seq, (*it).second)) {
			it = dump.erase(it);
		} else it++;
	}
	
	printf("after first filter: %lu sequences\n", dump.size());
	
	//второй фильтр - по очкам----------------------------------------------------
	for (auto it = dump.begin(); it != dump.end(); ) {
		if (filter2(p_count, score_matrix, index_arr, alphabet.length(), 
				seq, (*it).second)) {
			it = dump.erase(it);
		} else it++;
	}
	
	printf("after second filter: %lu sequences\n", dump.size());
	
	//Смит-Ватерман для "финалистов"
	std::vector<SWres> results;
	while(!dump.empty()) {
		results.push_back(SmithWaterman(score_matrix, index_arr, alphabet.length(), 
											penalty, seq, dump.begin()->second));
		
		//получить имяпоследовательности из базы
		sqlite3_stmt *stmt;
		const char *pzTest;
		std::string select = "SELECT description FROM Sequence WHERE ID = ?";
		
		rc = sqlite3_prepare_v2(db, select.c_str(), select.length(), 
														&stmt, &pzTest);
		if( rc == SQLITE_OK ){
			sqlite3_bind_int(stmt, 1, dump.begin()->first);
			sqlite3_step(stmt);
			results.back().name = (char*)sqlite3_column_text(stmt, 0);
			sqlite3_finalize(stmt);
		} else { 
			printf("Can't prepare sqlite3_prepare_v2 'SELECT description \
							FROM Sequence WHERE ID = %d'\n", dump.begin()->first
						); 
		}
		
		dump.erase(dump.begin());
	}
	
	//сортировка результатов
	std::sort(results.begin(), results.end());
	
	//вывод ответа
	int position = 0;
	for (auto it = results.begin(); it != results.end(); it++) {
		std::string res1 = "", res2 = "";
		GetAllign(it->way, it->i, it->j, it->seq1, it->seq2, res1, res2); 
		printf("\n#%d %s\n", ++position, it->name.c_str());
		printf(">score: %d\n", it->score);
		printf(">alignment:\n");
		printf("%s\n", res1.c_str());
		printf("%s\n", res2.c_str());
		for (int i = 0; i < res1.length(); i++) {
			char c = '!';
			if (res1[i] == ' ' || res2[i] == ' ' || res1[i] == res2[i]) c = ' ';
			printf("%c", c);
		}
		printf("\n");
		it->FreeSharedMem();
	}
	
	//free========================================================================
	delete[] score_matrix;
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