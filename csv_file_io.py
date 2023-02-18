"""
BSD 3-Clause License

Copyright (c) 2023, Özlem Salman
Copyright (c) 2023, Armağan Salman

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


#(
# Python system modules:
import csv
import typing as T
#)

#(
Bool = bool
Str = str
#)


def write_data_to_csv(file_path: Str \
		, data_rows: T.List[T.List[Str]] \
		, write_mode: Str = 'a' \
		, encoding: Str = 'utf-8' \
		, newline: Str = '') \
		-> Bool:
#(
	with open(file_path, write_mode, encoding=encoding \
			, newline=newline) as file_handle:
	#(
		csvwriter = csv.writer(file_handle)
		# Each row contains 0 or more elements (for example, a List).
		for row in data_rows:
		#(
			csvwriter.writerow(row)
		#)
	#)
	
	return True
#)


def read_data_from_csv(file_path: Str \
		, encoding: Str = 'utf-8' \
		, newline: Str = '') \
		-> T.List[T.List[Str]]:
#(
	with open(file_path, mode = 'r', encoding=encoding \
			, newline=newline) as file_handle:
	#(
		csvreader = csv.reader(file_handle)
		row_list = []
		
		# Each row contains 0 or more elements (for example, a List).
		for row in csvreader:
		#(
			row_list.append(row)
		#)
	#)
	
	return row_list
#)


def main(args):
#(
	data_rows = \
	[
		  ["Line 1"] \
		, ["Line 2"] \
		, ["Line 3"] \
		, ["Line 4"] \
	]
	
	file_path = "./CSV-test-0.11159914797446868.csv"
	
	success_bool = write_data_to_csv(file_path, data_rows, write_mode='w')
	
	if success_bool:
	#(
		print(f"Row count: {len(data_rows)}")
		print(f"Data written successfully to: {file_path}")
	#)
	
	read_data = read_data_from_csv(file_path)
	
	second_file_path = "./CSV-test-0.11159914797446868-2.csv"
	
	success_bool = write_data_to_csv(second_file_path, read_data, write_mode='w')
	
	if success_bool:
	#(
		print(f"Row count: {len(read_data)}")
		print(f"Data written successfully to: {second_file_path}")
	#)
	
	assert(len(read_data) == len(data_rows))
#)


if __name__ == "__main__":
#(
	args = {}
	
	main(args)
#)






