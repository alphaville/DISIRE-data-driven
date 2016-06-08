classdef DisireConnection < Singleton
    %DISIRECONNECTION establishes a connection to the DISIRE data server
    %from which one may obtain data for the walking beam furnace.
    
    properties(Access=private)
        myConnection;
    end
    
    properties(Constant,Access=private)
        WBF_RAW_DATA = 'WBFData';
        WBF_FINE_DATA = 'WBFDataShort';
    end
    
    methods(Access=private)
        function newObj = DisireConnection()
            host = '147.102.82.32:3306';
            user = 'monty';
            password = 'abbamammamia';
            dbName = 'WBF';
            jdbcString = sprintf('jdbc:mysql://%s/%s', host, dbName);
            jdbcDriver = 'com.mysql.jdbc.Driver';
            dbConn = database(dbName, user, password, jdbcDriver, jdbcString);
            newObj.myConnection = dbConn;
        end
    end
    
    methods(Static)
        function obj = instance()
            persistent uniqueInstance
            if isempty(uniqueInstance)
                obj = DisireConnection();
                uniqueInstance = obj;
            else
                obj = uniqueInstance;
            end
        end
    end
    
    methods % Public Access
        function idx = getRawData(obj)
            idx=[];
        end
        
        function idx = getFineData(obj, k1, k2, property)
            cmd_tmpl = 'select * from `%s` WHERE property="%s" AND `k` BETWEEN %d AND %d ORDER BY `k`';
            cmd = sprintf(cmd_tmpl, obj.WBF_FINE_DATA, property, k1, k2);
            idx = fetch(obj.myConnection, cmd);
        end                
        
        function delete(obj)
            close(obj.myConnection);
        end
    end
    
end