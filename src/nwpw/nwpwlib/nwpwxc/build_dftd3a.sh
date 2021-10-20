#!/usr/bin/env bash
do_exit(){
	echo ' '
	echo 'please install the patch command'
	echo ' '
	exit 1
    }
rm -f dftd3.f nwpwxc_vdw3a.F
export PATH=`pwd`:$PATH
if [[ ! -x "$(command -v patch)" ]]; then
    #try to download busybox for x86_64 linux
    if [[ $(uname -s) == "Linux" ]]  && [[ $(uname -m) == "x86_64" ]] ; then
	echo "patch command missing"
	echo "downloading busybox to use patch command"
	wget https://www.busybox.net/downloads/binaries/1.31.0-defconfig-multiarch-musl/busybox-x86_64 -O patch
	if [ "$?" != 0 ]; then
	    do_exit
	else
	    chmod +x patch
	fi
    else
	do_exit
    fi
fi
URL1="https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/"
URL2="https://web.archive.org/web/20210527062154if_/https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/"
declare -a urls=("$URL1" "$URL1" "$URL1" "$URL2" "$URL2")
TGZ=dftd3.tgz
# check gzip integrity
[ -f "$TGZ" ] && gunzip -t "$TGZ" > /dev/null && REUSE_TGZ=1
if [[ "${REUSE_TGZ}" == 1 ]]; then
    echo "using existing" "$TGZ"
else
    rm -f "$TGZ"
    echo "downloading"  "$TGZ"
    if [ "$(command -v curl)" ]; then
	COM1="curl -f -L"
	COM2="-o $TGZ"
    elif [ "$(command -v wget)" ]; then
	COM1="wget"
	COM2=" "
    else
	echo 'ERROR'
	echo 'please install either curl or wget to continue '
	exit 1
    fi
    echo chchc ${COM1} ${urls[$tries]}/"$TGZ" ${COM2}
    tries=1 ; until [ "$tries" -ge 5 ] ; do
    ${COM1} ${urls[$tries]}/"$TGZ" ${COM2} && break
    tries=$((tries+1)) ; echo attempt no.  $tries    ; sleep 5 ;  done
fi
if [[ ! -f "$TGZ" ]]; then
    echo "download failed"
    exit 1
fi
tar xzf dftd3.tgz dftd3.f
mv dftd3.f nwpwxc_vdw3a.F
patch -p0 < nwpwxc_vdw3a.patch

