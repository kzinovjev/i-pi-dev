#!/bin/bash
DRIVER=/ssd/i-pi-ghts/bin/i-pi-driver
IPI=/ssd/i-pi-ghts/bin/i-pi

i=$((OMPI_COMM_WORLD_RANK+1))
nproc=$OMPI_COMM_WORLD_SIZE
nbeads=32
temp=300
port=$((21142+i))

if [ $i -le $(($nproc/2)) ];
then
    xyz="react"
else
    xyz="prod"
fi

rm -rf $i
mkdir $i
cd $i

cp ../../eckart.xml xml
sed -i "s/__XYZ__/$xyz/g" xml
sed -i "s/__PORT__/$port/g" xml
sed -i "s/__NBEADS__/$nbeads/g" xml
sed -i "s/__SEED__/$RANDOM/g" xml
sed -i "s/__TEMP__/$temp/g" xml

cp ../../ghts.json .
sed -i "s/__WALKER__/$i/g" ghts.json

$IPI xml &
echo $! >> ../../pids
while [ ! -e /tmp/ipi_localhost$port ]; do
    sleep 1
done
$DRIVER -u -h localhost$port -m eckart -p $port -o -0.5,2,0.3,0.551152161 > eckart.out &
echo $! >> ../../pids
wait
echo "finished walker $i"
