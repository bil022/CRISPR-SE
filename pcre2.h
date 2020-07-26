#ifndef PCRE2_H
#define PCRE2_H

#define PCRE2_CODE_UNIT_WIDTH 8

#include <stdio.h>
#include <iostream>
#include <string>
#include <pcre2.h>
#include <assert.h>

using namespace std;

class Pcre2 {
private:
	pcre2_code *re;
	PCRE2_SPTR pattern;     /* PCRE2_SPTR is a pointer to unsigned code units of */
	PCRE2_SPTR subject;     /* the appropriate width (in this case, 8 bits). */
	
	int crlf_is_newline, errornumber, rc, utf8;
	uint32_t option_bits, newline;
	PCRE2_SIZE erroroffset, *ovector;
	size_t subject_length;
	pcre2_match_data *match_data;
public:
	Pcre2() {
		match_data=NULL;
		re=NULL;
	}
	~Pcre2() {
		if (match_data) {
			pcre2_match_data_free(match_data);
			match_data=NULL;
		}
		if (re) {
			pcre2_code_free(re);
			re=NULL;
		}
	}
	virtual void found(int offset) = 0;
	int regex(const char* pam, string& seq) {
		pattern = (PCRE2_SPTR)pam;
		subject = (PCRE2_SPTR)seq.c_str();
		subject_length = seq.length();
		
		re = pcre2_compile(
						   pattern,               // the pattern
						   PCRE2_ZERO_TERMINATED, // indicates pattern is zero-terminated
						   0,                     // default options
						   &errornumber,          // for error number
						   &erroroffset,          // for error offset
						   NULL);                 // use default compile context
		
		if (re == NULL) {
			PCRE2_UCHAR buffer[256];
			pcre2_get_error_message(errornumber, buffer, sizeof(buffer));
			cerr << "PCRE2 compilation failed at offset " << erroroffset << ": " << buffer << endl;
			return 1;
		}
		
		match_data = pcre2_match_data_create_from_pattern(re, NULL);
		
		rc = pcre2_match(
						 re,                   /* the compiled pattern */
						 subject,              /* the subject string */
						 subject_length,       /* the length of the subject */
						 0,                    /* start at offset 0 in the subject */
						 0,                    /* default options */
						 match_data,           /* block for storing the result */
						 NULL);                /* use default match context */
		
		if (rc < 0) {
			switch(rc) {
				case PCRE2_ERROR_NOMATCH: cerr << "No match\n"; break;
				default: cerr << "Matching error " << rc << endl; break;
			}
			return 1;
		}
		
		ovector = pcre2_get_ovector_pointer(match_data);
		found(ovector[0]);
		
		if (rc == 0) {
			cerr << "ovector was not big enough for all the captured substrings\n";
		}
		
		if (ovector[0] > ovector[1]) {
			cerr << "Run abandoned\n";
			return 1;
		}
		
		(void)pcre2_pattern_info(re, PCRE2_INFO_ALLOPTIONS, &option_bits);
		utf8 = (option_bits & PCRE2_UTF) != 0;
		
		(void)pcre2_pattern_info(re, PCRE2_INFO_NEWLINE, &newline);
		crlf_is_newline = newline == PCRE2_NEWLINE_ANY ||
		newline == PCRE2_NEWLINE_CRLF ||
		newline == PCRE2_NEWLINE_ANYCRLF;
		
		for (;;)
		{
			PCRE2_SIZE start_offset = ovector[0]+1;
			assert(ovector[0] != ovector[1]);
			
			PCRE2_SIZE startchar = pcre2_get_startchar(match_data);
			if (start_offset <= startchar)
			{
				if (startchar >= subject_length) break;   /* Reached end of subject.   */
				start_offset = startchar + 1;             /* Advance by one character. */
				if (utf8) {                                 /* If UTF-8, it may be more  */
					for (; start_offset < subject_length; start_offset++)
						if ((subject[start_offset] & 0xc0) != 0x80) break;
				}
			}
			
			rc = pcre2_match(
							 re,                   /* the compiled pattern */
							 subject,              /* the subject string */
							 subject_length,       /* the length of the subject */
							 start_offset,         /* starting offset in the subject */
							 0,              /* options */
							 match_data,           /* block for storing the result */
							 NULL);                /* use default match context */
			
			if (rc == PCRE2_ERROR_NOMATCH) {
				break;                    /* All matches found */
			}
			
			if (rc < 0) {
				cerr << "Matching error " << rc << endl;
				return 1;
			}
			
			found(ovector[0]);
			
			if (rc == 0) {
				cerr << "ovector was not big enough for all the captured substrings\n";
			}
			
			if (ovector[0] > ovector[1]) {
				cerr << "Run abandoned\n";
				return 1;
			}
		}
		
		return 0;
	}
};

#endif
