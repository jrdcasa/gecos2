import sqlite3
import sys

class DBjobs(object):

    """Class to represent the database of jobs

    This class is used to keep track the QM jobs in the server.
    The database has a table called ``qm_jobs`` with the following columns:

        * ``id_job`` (integer): Primary key of the table
        * ``name_job`` (string): Name of the input file in the QM package.\
        The pattern is **<name_of_the_qm_package>_<identification>_<0,1,2>_<1,2>**\
        Example: nwchem_00024_0_1 (for segment 1) or gauss_00374_2_1 (for pair 2,1)
        * ``status_job`` (string): Status of the job in the server
        * ``pid_slurm`` (integer): Job PID from slurm queue management
        * ``energy`` (float): Electronic energy from the output
        * ``cog`` (list[3]): Center of geometry

    **Attributes**

    Attributes:
        _con (Connection): Connecting to the SQLite database (sqlite3.Connection)
        _cursor (Cursor): Connection object to execute SQLite queries from Python. (sqlite3.Cursor)
        _dbname (str): Path to the database file
        _logger (Logger): Logger to send messages

    **Methods**

    """

    # ***********************************************************************************
    def __init__(self, dbname, logger=None):

        """Constructor of a DBjobs object

        Args:
            dbname (str): Path to the database
            logger (logger): Logger to send messages

        """

        self._dbname = dbname
        self._logger = logger

        # Create connection to the database
        try:
            self._con = sqlite3.connect(dbname)
            self._cursor = self._con.cursor()
        except sqlite3.Error:
            print(sqlite3.Error)

    # ************************************************************************************
    def check_for_table(self, tablename):

        """Checks if table ``tablename`` exists

        Args:
            tablename (string): Name of the table in the database

        Returns:
            ``True`` if the tablename exists otherwise ``False``

        """

        try:
            is_table = self._db_schema().index(tablename)
            return True
        except ValueError:
            return False

    # ***********************************************************************************
    def close_connection(self):

        """Close the database connection

        """

        self._con.close()

    # ***********************************************************************************
    def commit_db(self):

        """Commits the database
        """

        self._con.commit()

    # ***********************************************************************************
    def create_table_qmjobs(self):

        """ Create table ``qmjobs`` if it does not exist in the database

        Returns:
            ``True`` if the table is created otherwise ``False``

        """

        try:
            self._cursor.execute("CREATE TABLE qm_jobs(id_job integer PRIMARY KEY, \
                                                          name_job text, \
                                                          status_job text, \
                                                          pid_slurm integer, \
                                                          energy real, \
                                                          cog text, \
                                                          cpu_time real)")
            success = True
        except sqlite3.Error:
            if self._logger is None:
                #print("The table qmjobs in database {} already exists".format(self._dbname))
                success = False
            else:
                #self._logger.info("\tThe table qmjobs in database {} already exists".format(self._dbname))
                success = False
            pass

        self._con.commit()
        return success

    # ***********************************************************************************
    def insert_data(self, sql):

        """Insert data in the datababase using a valid SQL sentence

        If data does not exist, then insert into the database otherwise not insertion is done

        Using a sql INSERT sentence

        Args:
            sql (string): A valid sql sentence

        Returns:
            ``True`` if data are inserted otherwise ``False``

        Example:
            sql_insert = "INSERT INTO qm_jobs VALUES ({0:d}, '{1:s}', 'CREATED', 'NULL', 'NULL', '{2:s}')".\
            format(self._db_index, inputname, str(cog_list)[1:-1])
            db.insert_data(sql_insert):
        """

        try:
            self._cursor.execute(sql)
            return True
        except sqlite3.Error:
            #print(sys.exc_info())
            return False

    # ***********************************************************************************
    def query_data(self, sql):

        try:
            return self._cursor.execute(sql)
        except sqlite3.Error as e:
            print(e)
            return None

    # ************************************************************************************
    def number_of_rows(self, tablename):

        """Returns the number of rows of the table

        Args:
            tablename (str): Name of the table

        Returns:
            (int) Number of rows

        """

        sql = '''SELECT COUNT(*) from {0:}'''.format(tablename)
        self._cursor.execute(sql)

        nrows = self._cursor.fetchone()

        return nrows[0]

    # ***********************************************************************************
    def remove_row(self, tablename, fieldname, id_job):

        """Remove a row

        Args:
            tablename (string): Name of hte table
            fieldname (string): Name of the column to delete
            id_job (any): Value for the fieldname

        Returns:
            ``True`` if the row is remove otherwise ``False``

        """

        try:
            sql = 'DELETE FROM {0:} WHERE {1:}=?'.format(tablename, fieldname)
            self._cursor.execute(sql, (id_job,))
            self.commit_db()
            return True
        except sqlite3.Error:
            #print(sys.exc_info())
            return False

    # ***********************************************************************************
    def remove_table(self, tablename):

        """Remove table in the database

        Args:
            tablename (str):

        Returns:

        """

        sql = """DROP TABLE {0:}""".format(tablename)
        self._cursor.execute(sql)

    # ***********************************************************************************
    def sort_by_column(self, tablename, columnname):

        """Sort columns

        Args:
            tablename (str): Name of the table
            columnname (str): Column name in the table

        Returns:
            (tuple), list with the ordered rows, maximum index

        """

        sql = "SELECT * FROM {} ORDER BY {}".format(tablename, columnname)
        rows = self._cursor.execute(sql)
        rows_list = rows.fetchall()
        index_max = rows_list[0][0]
        return rows_list, index_max

    # ***********************************************************************************
    def update_allcolumns_tosamevalue(self, table, d):

        """Update all rows of a certain column in the dictionary

        Args:
            table (str): Name of the table
            d (dict): Dictionary {<name_of_column>: <value for the whole column>}

        Returns:
            True if the operation is success.

        """

        try:
            for key, value in d.items():
                sql = """UPDATE {} SET {} = '{}'""".format(table, key, d[key])
                self._cursor.execute(sql)
            return True
        except sqlite3.Error:
            #print(sys.exc_info())
            return False

    # ***********************************************************************************
    def update_data_row(self, table, field, field_value, where_field, where_value):

        """Update date by row using where

        Args:
            table (str):
            field (str):
            field_value (any):
            where_field (str):
            where_value (any):

        Returns:

        """

        sql = """UPDATE {0:} SET {1:} = '{2:}' WHERE {3:} = '{4:}'""".\
            format(table, field, field_value, where_field, where_value)
        try:
            self._cursor.execute(sql)
        except:
            print("CJJJJ")
            exit()

    # ************************************************************************************
    def _db_schema(self):

        """Show the name of the tables in the database

        """


        sql = "SELECT name FROM sqlite_master WHERE type='table'"
        self._cursor.execute(sql)
        list_of_tuples = self._cursor.fetchall()
        l = list(map(lambda x: x[0], list_of_tuples))
        return l

    # ************************************************************************************
    def account_qm(self, tablename):

        """
        Account for the total number of QM and the completed QM jobs in the database
        The format of the lists are:

        total_qm_pairs    = [2,2,2,2]
        complete_qm_pairs = [1,1,1,1] --> [segment_1_1, segment_1_2, segment_2_1, segment_2_2]

        Args:
            tablename (str): Name of the table to find

        Returns:
            A tuple of two lists. The total_qm_pairs and complete_qm_pairs

        """

        # Pair 1-0, 2-0, 1-1, 1-2, 2-1, 2-2
        total_qm_pairs = []
        complete_qm_pairs = []
        running_qm_pairs = []
        inserver_qm_pairs = []
        pending_qm_pairs = []
        error_qm_pairs = []

        for iseg in range(1,3):
            sql = 'SELECT * FROM {0:} WHERE name_job LIKE \"%{1:1d}_{2:1d}\"'. \
                format(tablename, 0, iseg)
            self._cursor.execute(sql)
            rows_total = self._cursor.fetchall()
            total_qm_pairs.append(len(rows_total))

            sql_completed = 'SELECT * FROM {0:} WHERE name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "COMPLETED" '. \
                format(tablename, 0, iseg)
            self._cursor.execute(sql_completed)
            rows = self._cursor.fetchall()
            complete_qm_pairs.append(len(rows))

            sql_error = 'SELECT * FROM {0:} WHERE ' \
                        '(name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "RUNNING") OR ' \
                        '(name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "PENDING") OR ' \
                        '(name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "INSERVER") OR ' \
                        '(name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "COMPLETED")  '. \
                format(tablename, 0, iseg)
            self._cursor.execute(sql_error)
            rows = self._cursor.fetchall()
            error_qm_pairs.append(len(rows_total) - len(rows))

            sql_running = 'SELECT * FROM {0:} WHERE name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "RUNNING" '. \
                format(tablename, 0, iseg)
            self._cursor.execute(sql_running)
            rows = self._cursor.fetchall()
            running_qm_pairs.append(len(rows))

            sql_pending = 'SELECT * FROM {0:} WHERE name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "PENDING" '. \
                format(tablename, 0, iseg)
            self._cursor.execute(sql_pending)
            rows = self._cursor.fetchall()
            pending_qm_pairs.append(len(rows))

            sql_inserver = 'SELECT * FROM {0:} WHERE name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "INSERVER" '. \
                format(tablename, 0, iseg)
            self._cursor.execute(sql_inserver)
            rows = self._cursor.fetchall()
            inserver_qm_pairs.append(len(rows))

        for iseg in range(1,3):
            for jseg in range(1,3):
                sql = 'SELECT * FROM {0:} WHERE name_job LIKE \"%{1:1d}_{2:1d}\"'.\
                    format(tablename, iseg, jseg)
                self._cursor.execute(sql)
                rows_total = self._cursor.fetchall()
                total_qm_pairs.append(len(rows_total))

                sql_completed = 'SELECT * FROM {0:} WHERE name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "COMPLETED" '.\
                    format(tablename, iseg, jseg)
                self._cursor.execute(sql_completed)
                rows = self._cursor.fetchall()
                complete_qm_pairs.append(len(rows))

                sql_error = 'SELECT * FROM {0:} WHERE ' \
                            '(name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "RUNNING") OR ' \
                            '(name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "PENDING") OR ' \
                            '(name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "INSERVER") OR '\
                            '(name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "COMPLETED")  '. \
                            format(tablename, iseg, jseg)
                self._cursor.execute(sql_error)
                rows = self._cursor.fetchall()
                error_qm_pairs.append(len(rows_total)-len(rows))

                sql_running = 'SELECT * FROM {0:} WHERE name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "RUNNING" '. \
                    format(tablename, iseg, jseg)
                self._cursor.execute(sql_running)
                rows = self._cursor.fetchall()
                running_qm_pairs.append(len(rows))

                sql_pending = 'SELECT * FROM {0:} WHERE name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "PENDING" '. \
                    format(tablename, iseg, jseg)
                self._cursor.execute(sql_pending)
                rows = self._cursor.fetchall()
                pending_qm_pairs.append(len(rows))

                sql_inserver = 'SELECT * FROM {0:} WHERE name_job LIKE \"%{1:1d}_{2:1d}\" AND status_job == "INSERVER" '. \
                    format(tablename, iseg, jseg)
                self._cursor.execute(sql_inserver)
                rows = self._cursor.fetchall()
                inserver_qm_pairs.append(len(rows))

        return total_qm_pairs, inserver_qm_pairs,  pending_qm_pairs, complete_qm_pairs, error_qm_pairs, running_qm_pairs

    # ************************************************************************************
    def add_column(self, tablename, columname, columntype):

        """

        Args:
            tablename:
            columname:
            columntype:

        Returns:

        """

        try:
            sql = "ALTER TABLE {} ADD COLUMN {} {}".format(tablename, columname, columntype)
            self._cursor.execute(sql)
        except sqlite3.OperationalError:
            pass
