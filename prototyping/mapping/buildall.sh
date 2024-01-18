cd tinyxml; echo -e "\n###building tinyxml"; make -j install; cd ..;
cd ann_1.1.2; echo -e "\n###building ann_1.1.2"; make -j install; cd ..;
cd GiDParser; echo -e "\n###building GiDParser"; make -j install; cd ..;
cd mortarMapper; echo -e "\n###building mortarMapper"; make -j install; cd ..;
cd dataCreator; echo -e "\n###building dataCreator"; make -j install; cd ..;
cd testMapper; echo -e "\n###building testMapper"; make -j