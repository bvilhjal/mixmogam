import env
import sys

user_name = env.env['db_user']
user_passwd = env.env['db_passwd']

def connect_to_papaya(db="stock_250k"):
	return connect_to_db("papaya.usc.edu", db, user_name, user_passwd)

def connect_to_gmi_ara_devel_be(db="stock_250k"):
	return connect_to_db("gmi-ara-devel-be.gmi.oeaw.ac.at", db, user_name, user_passwd)

def connect_to_default_lookup(db="stock_250k"):
	import env
	return connect_to_db(env.env['default_lookup_db'], db, user_name, user_passwd)

def connect_to_default_insert(db="stock_250k"):
	import env
	return connect_to_db(env.env['default_insert_db'], db, user_name, user_passwd)

def connect_to_arabidopsis(db='stock_250k', port=3309):
	return connect_to_db('arabidopsis.gmi.oeaw.ac.at', db, user_name, user_passwd, port=port)


def connect_to_db(host, db, user=user_name, passwd=user_passwd, port=3306):
	"""
	Returns a MySQLdb connection (needs to be closed).
	"""
	import MySQLdb
	#Load cand. gene list.	
	print "Connecting to db, host=%s, db=%s " % (host, db)
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db="stock_250k", port=port)
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	return conn




def get_tg_ecotypes(ecotypes, user=None, passwd=None, host="papaya.usc.edu"):

	import MySQLdb
	print "Connecting to db, host=" + host
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host=host, user=user, passwd=passwd, db="at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()

	url = salc.engine.url.URL("mysql", username=user_name, password="bamboo123", host="papaya.usc.edu", database="stock")
	e = salc.create_engine(url) #Creating an engine object
	#c = e.connect()
	meta = salc.MetaData()
	meta.bind = e
	einfo_table = salc.Table('ecotype_info', meta, autoload=True)
	from sqlalchemy.sql import select
	res = select([einfo_table])
	print res

if __name__ == "__main__":
	get_tg_ecotypes(1)
