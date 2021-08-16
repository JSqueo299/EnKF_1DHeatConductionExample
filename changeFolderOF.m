function changeFolderOF(varname,tFolder,caseFolder_OF)
    
    if isreal(tFolder) && tFolder >= 0
        str1 = sprintf('%s/%.15g',caseFolder_OF,tFolder);
        cd(str1);
    else
        error('Cannot change to specified time folder.\n Time Folder is not a real number or is a negative value');
    end

end

