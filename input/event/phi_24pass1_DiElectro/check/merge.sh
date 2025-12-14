echo ${list}
for cluster in `seq 1 4`; do
  for period in af ag aj al am an ao; do
    echo ${cluster}/${period}
    #grep -Ff list/cluster_1_985_inRow.list list/file_local_24pass1 | grep -f list/period/af
    list=`grep -Ff list/cluster_${cluster}*_inRow.list list/file_local_24pass1 | grep -f list/period/${period}`
    list=`echo "${list}" | sort | uniq`
    if [ -z "${list}" ]; then
      echo "No files for cluster ${cluster} period ${period}"
      continue
    fi
    list=`echo "${list}" | tr '\n' ' '`
    echo hadd -f merge/${cluster}_${period}.root ${list}
    hadd -f merge/${cluster}_${period}.root ${list}
  done
done
