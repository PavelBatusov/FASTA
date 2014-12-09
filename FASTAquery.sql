PRAGMA foreign_keys = ON;
PRAGMA foreign_keys;

DROP TABLE if exists SubSequence;
DROP TABLE if exists Sequence;

CREATE TABLE Sequence (
	ID INTEGER NULL PRIMARY KEY AUTOINCREMENT,
	description TEXT,
	string TEXT
);

CREATE TABLE SubSequence (
	ID INTEGER NULL PRIMARY KEY AUTOINCREMENT,
	string TEXT,
	position INTEGER,
	baseString INTEGER,
	FOREIGN KEY(baseString) REFERENCES Sequence(ID) ON DELETE CASCADE
);

$ 