/**
 * @file ErrorLogger.cpp
 *
 * @brief A simple threadsafe singleton error logger that stores and counts 
 * errors.
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <mutex>

using namespace std;

const string EMPTY_TYPE = "";

struct ErrorData {

	ErrorData(const string &msg, const string &type) : msg(msg), type(type) {}

	const string msg;
	const string type;

};

class ErrorLogger {

public:

	ErrorLogger   (      ErrorLogger &other) = delete;
	void operator=(const ErrorLogger &other) = delete;

	void clear            (                                                   );
	void logError         (const string  &msg, const string &type = EMPTY_TYPE);
	void setOutputStream  (      ostream &out                                 );

	size_t countErrors(const string &type = EMPTY_TYPE) const;

	static ErrorLogger &getInstance();

private:

	ErrorLogger();

	vector<ErrorData> errors;

	ostream *errorStream;

	mutable mutex errorLock;

};

ErrorLogger::ErrorLogger() : errorStream(&cerr) {}

void ErrorLogger::clear() {

	errorLock.lock();
	errors.clear();
	errorLock.unlock();

}

void ErrorLogger::logError(
	const string &msg, 
	const string &type = EMPTY_TYPE
) {

	errorLock.lock();
	*errorStream << msg << endl;
	errors.emplace_back(msg, type);
	errorLock.unlock();

}

void ErrorLogger::setOutputStream(ostream &out) {

	errorLock.lock();
	errorStream = &out;
	errorLock.unlock();

}

size_t ErrorLogger::countErrors(const string &type = EMPTY_TYPE) const {

	size_t val;

	errorLock.lock();

	if(type == EMPTY_TYPE) {

		val = errors.size();

	} else {

		val = count_if(
			errors.cbegin(),
			errors.cend(),
			[type](const ErrorData &e) {
				return e.type == type;
			}
		);

	}

	errorLock.unlock();

	return val;

}

ErrorLogger &ErrorLogger::getInstance() {

	static ErrorLogger instance;

	return instance;

}