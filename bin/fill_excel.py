#!/usr/bin/env python3
import sys
from openpyxl import load_workbook
#from openpyxl.utils.dataframe import dataframe_to_rows

def main(argv):#argv):
   row_start = 4
   wb = load_workbook(argv[-1])
   cov_list = argv[1:-2]
   ws = wb.active
   for ix,seq_name in enumerate(cov_list[1::2]):
      ws["A"+str(row_start)] = "submitter_name" 
      ws["C"+str(row_start)] = seq_name 
      ws["D"+str(row_start)] = "betacoronavirus"
      ws["E"+str(row_start)] = "original"
      ws["I"+str(row_start)] = "Human"
      ws["R"+str(row_start)] = "Nanopore MinION"
      ws["S"+str(row_start)] = "ARTIC-nCoV-bioinformaticsSOP-v1.1.0"
      ws["T"+str(row_start)] = cov_list[ix*2]
      ws["U"+str(row_start)] = "Molecular Medicine Laboratory, University of Magallanes"
      ws["V"+str(row_start)] = "Av. los Flamencos 01364"
      ws["X"+str(row_start)] = "Centro Asistencial Docente y de Investigacion, Universidad de Magallanes"
      ws["Y"+str(row_start)] = "Av. los Flamencos 01364"
      ws["AA"+str(row_start)] = "Jorge González, Jacqueline Aldridge, Diego Alvarez, Ines Cid, Constanza Ceron, Roberto Uribe-Paredes, Marcelo Navarrete"
      row_start = row_start+1
   

   wb.save("Submitted.xlsx")
if __name__ == "__main__":
    main(sys.argv[1:])
