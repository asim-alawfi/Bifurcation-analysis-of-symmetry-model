#!/usr/bin/python
#
# $Id: c_insert.py 91 2015-01-11 17:18:05Z jansieber $
#
import os, sys, re
import subprocess
usage='''\
Usage: python {0:s} rootdir version
Includes version and (c) info line in all m and html files. 
If Id line is present, then Id line is replaced with (c) line using argument version
'''

#matching patterns
filetypes=['.*[.]m']

# convert matched regexp to sortable date integer
getdate=lambda x, order:int(x[order[2]])*10000+int(x[order[1]])*100+int(x[order[0]])

# main function performing all steps
def main():
    nargs=len(sys.argv)
    if nargs<3:
        print(usage.format(sys.argv[0]))
        print('nargs:', nargs)
        print('argv:', sys.argv)
        return -1
    version=sys.argv[2]
    rootdir=sys.argv[1]
    files=filelist(rootdir, filetypes)
    for k in range(len(files)):
        f=open(files[k], 'r')
        text=f.readlines()
        f.close()
        cmd_str=["git","log","-1","--format=%cd (%H)","--date=short",files[k]]
        lastver=subprocess.Popen(cmd_str,stdout=subprocess.PIPE).stdout.read()
        lastver=lastver.decode("utf-8")
        insert_str="% (c) DDE-Biftool v"+version+" "+lastver+"\n"
        #print(insert_str)
        firstcomment=find_first(r"^\w*%",text)
        firstId=find_first(r"^.*[$]Id",text)
        first_c=find_first(r"^.*%.*[(]c[)]",text)
        insert_num=-1
        insertframe=False
        if firstId>=0:
            text.pop(firstId)
            insert_num=firstId
        if first_c>=0:
                text.pop(first_c)
                if insert_num<0:
                    insert_num=first_c
        if firstcomment>=0:
            if insert_num<0:
                insert_num=firstcomment
                insertframe=True
        if insert_num<0:
            insert_num=0
        if insertframe:
            text.insert(insert_num,"%\n")
        text.insert(insert_num,insert_str)
        if insertframe:
            text.insert(firstcomment,"%\n")
        f=open(files[k], 'w')
        f.writelines(text)
        f.close()
        #print(files[k], insert_num)
    return 0
# find first occurence of regex in lst
def find_first(rg,lst):
    bool=[a is not None for a in (re.match(rg,c) for c in lst)]
    first=-1
    for i in range(len(bool)):
        if bool[i]:
            first=i;
            break
    return first
# recursive list of files with names matching fpatterns in path
def filelist(path, fpatterns):
    '''
    returns a list of all files in path (recursively) with filename pattern fpattern and 
    containing a line matching regex
    list elements have format [filename_with_path, number_of_first_matching_line]
    '''
    fp_re=[re.compile(k) for k in fpatterns]
    res = []
    for root, dirs, fnames in os.walk(path):
        for fname in fnames:
            for fp in fp_re:
                if fp.match(fname) is None:
                    continue
                res.append(os.path.join(root,fname));
    return res

# create (c) string with version (as string), commit number (if not None) and date
def c_create(version,  date, commit=None):
    year=date/10000
    month=(date-year*10000)/100
    day=date-year*10000-month*100
    if commit:
        s='(c) DDE-BIFTOOL v. {0:s}({1:d}), {2:02d}/{3:02d}/{4:04d}'.format(version, commit, day, month, year)
    else:
        s='(c) DDE-BIFTOOL v. {0:s}, {2:02d}/{3:02d}/{4:04d}'.format(version, commit, day, month, year)
    return s

# scipt body
if __name__ == '__main__':
    sys.exit(main())
