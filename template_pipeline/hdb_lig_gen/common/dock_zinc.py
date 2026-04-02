#!/usr/bin/env python
'''connects to ZINC database for smiles, vendors, etc.

Michael Mysinger 200801 Created
Ryan Coleman 201204 Edited down to just SQL queries. See getposes.py for more.
'''

import os
import sys
import MySQLdb

REMARK_PREFIX = "#    "

# MySQL server
MYSQL_HOST = 'zincdb6'
MYSQL_USER = 'lab'
MYSQL_DATABASE = 'zinc8'

VENDOR_SQL = "select name, supplier_code from catalog_item join catalog " + \
    "on catalog_item.cat_id_fk=catalog.cat_id where sub_id_fk=%s"
SMILES_SQL = "select smiles from substance where sub_id=%s"

def init_zinc():
  db = MySQLdb.connect(host=MYSQL_HOST, user=MYSQL_USER, db=MYSQL_DATABASE)
  cursor = db.cursor()
  return db, cursor

def close_zinc(db, cursor):
  cursor.close()
  db.close()

def get_zinc_data(cursor, zid):
  retLines = []
  if zid[:4].isalpha():
    # ID format is ZINC12345678, trimming to 12345678 for lookup.
    zid = zid[4:]
  elif zid[:1].isalpha():
    # ID format is C12345678, trimming to 12345678 for lookup.
    zid = zid[1:]
  cursor.execute(SMILES_SQL, (zid))
  smiles = cursor.fetchone()
  if smiles:
    retLines.append(REMARK_PREFIX + "Substance Smiles   : %s\n" % smiles)
  cursor.execute(VENDOR_SQL, (zid))
  vendors = cursor.fetchall()
  for vendor in vendors:
    retLines.append(REMARK_PREFIX + "Vendor             : %s %s\n" % vendor)
  return retLines

if __name__ == "__main__":
  pass  # no way to run from commandline
