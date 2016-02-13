echo "init"
git init
find ./ -name '*.py' -or -name '*.pl' -or -name '*.R' -or -name '*.sh' | xargs git add
git commit -m "RIL figures"
git remote add origin https://github.com/JinfengChen/RILs_Preparation.git
git push -u origin master

echo "update"
find ./ -name '*.py' -or -name '*.pl' -or -name '*.R' -or -name '*.sh' | xargs git add
git commit -m "RIL figures 20160212"
git remote add update https://github.com/JinfengChen/RILs_Preparation.git
git push

