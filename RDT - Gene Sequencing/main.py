
# AUTHOR - RAJARSHI RAY

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from Bio import Entrez


print("---------- WELCOME TO DNA SEQUENCE TESTER AND DATABASE FOR STORING GENETIC INFORMATION ---------- \n")
print("DEVELOPED BY UEM KOLKATA BIOTECH DEPARTMENT BT-2 : RAJARSHI RAY")
print("----- KEY MAPPING ARE AS FOLLOWS ----- \n")
print(" \n A --> DNA ANALYSIS \n B --> CODON TABLE \n C --> CUSTOM DATABASE INPUT \n D --> ACCESS FASTA / GENBANK FILES IN SYSTEM \n E --> G-C CONTENT CALCULATOR \n F --> AMINO ACID MAP \n G --> ACCESS THE NCBI DATABASE \n Q --> EXIT THE APP \n")


# FUNCTION FOR TRANSCRIPTION AND TRANSLATION

def DNA_Analysis():
    inpSeq = input("PLEASE ENTER YOUR DNA SEQUENCE ---> ")
    mySeq = Seq(inpSeq)
    mySeq_complement = mySeq.complement()
    mySeq_mrna = mySeq.transcribe()
    proteinsSeq = mySeq_mrna.translate(to_stop=True)
    print(f"YOUR DNA SEQUENCE IS : {mySeq}")
    print(f"YOUR COMPLEMENT DNA SEQUENCE IS : {mySeq_complement}")
    print(f"YOUR mRNA SEQUENCE IS : {mySeq_mrna}")
    print(f"YOUR Protein SEQUENCE IS : {proteinsSeq}")

# FUNCTION FOR CODON TABLES


def Codon_Table():
    CODON_TABLE_SELECT = input(
        "WHICH CODON TABLE DO YOU WANT TO SEE ? [S : STANDARD || M : MITOCHONDRIAL] : ")
    if CODON_TABLE_SELECT == "S":
        getCodonTable_STANDARD = CodonTable.unambiguous_dna_by_name["Standard"]
        print(
            f"STANDARD CODON TABLE WITH RESPECTIVE AMINO ACIDS : \n {getCodonTable_STANDARD}")
    elif CODON_TABLE_SELECT == "M":
        getCodonTable_MITOCHONDRIAL = CodonTable.unambiguous_dna_by_id[2]
        print(
            f"MITOCHONDRIAL CODON TABLE WITH RESPECTIVE AMINO ACIDS : \n {getCodonTable_MITOCHONDRIAL}")

# FUNCTION FOR USER CUSTOM DATABASE


def Custom_Data():
    seqUser_INP = input("WHAT IS THE DNA SEQUENCE YOU HAVE RECORDED ? : ")
    seqUser = Seq(seqUser_INP)

    seqName = input("WHAT IS THE NAME OF YOUR RECORDED SUBJECT ? : ")
    seqDesc = input("ANY DESCRIPTION YOU WANT TO ADD ? : ")
    seqID = input("ID OF YOUR RECORD ? : ")
    seqFeat = input("FEATURES OF YOUR RECORD : ")
    seqFeat_FINAL = seqFeat.split()

    seqRECORD = SeqRecord(seqUser, name=seqName,
                          description=seqDesc, id=seqID, features=seqFeat_FINAL)

    print(f"YOUR FOLLOWING RECORD IS : \n {seqRECORD}")

    with open("data.txt", "a") as data:
        data.write(
            "NEW DATA VALUES ADDED : \n" "----------------- \n" f"DNA SEQUENCE IS : {seqUser} \n" f"SUBJECT NAME : {seqName} \n" f"ID OF SUBJECT : {seqID} \n" f"FEATURES : {seqFeat} \n" "----------------- \n" "END OF CURRENT DATABASE \n")

# FUNCTION FOR ACCESSING FASTA / GENBANK FILES


def Access_Files():
    print("YOU NEED TO STORE YOUR FASTA & GENBANK FILES IN LOCAL PARENT DIRECTORY : \n")
    fileName = input("NAME OF YOUR FILE : ")
    fileType = input("SPECIFY YOUR FILE TYPE [FASTA / GENBANK] : ")

    accessFile = SeqIO.read(fileName, fileType)
    print(f"FILE OPENED : \n {accessFile}")

# CHECK G-C CONTENT OF DNA SAMPLE


def GC_Content():
    inpSeq = input("PLEASE ENTER YOUR DNA SEQUENCE ---> ")
    mySeq = Seq(inpSeq)
    Adenine_Count = mySeq.count("A")
    Guanine_Count = mySeq.count("G")
    Thymine_Count = mySeq.count("T")
    Cytosine_Count = mySeq.count("C")

    GC_CONTENT = ((Guanine_Count + Cytosine_Count) / (Adenine_Count +
                  Thymine_Count+Guanine_Count+Cytosine_Count)) * 100

    print(f"THE G-C COUNT FOR THE DNA SEQUENCE : {mySeq} is = {GC_CONTENT} %")

# AMINO ACID MAPS


def Amino_Acid_Map():
    Amino_Acids = {
        "M": "METHIONINE",
        "A": "ALANINE",
        "R": "ARGININE",
        "N": "ASPARAGINE",
        "D": "ASPARTIC ACID",
        "C": "CYSTEINE",
        "Q": "GLUTAMINE",
        "E": "GLUTAMIC ACID",
        "G": "GLYCINE",
        "H": "HISTIDINE",
        "I": "ISOLEUCINE",
        "L": "LEUCINE",
        "K": "LYSINE",
        "F": "PHENLYALANINE",
        "P": "PROLINE",
        "O": "PYROLYSINE",
        "S": "SERINE",
        "U": "SELENOCYSTEINE",
        "T": "THREONINE",
        "W": "TRYPTOPHAN",
        "Y": "TYROSINE",
        "V": "VALINE",
        "B": "ASPARTIC ACID / ASPARAGINE",
        "Z": "GLUTAMIC ACID / GLUTAMINE",
    }
    [print(key, ":", value) for key, value in Amino_Acids.items()]

# ACCESS THE NCBI DATABASE AND RETRIEVE / DOWNLOAD DATA


def NCBI_Entrez():
    print("THIS IS JUST A TEST VERSION --- CAN ENCOUNTER ERRORS \n")
    print("--------------------------------------------------------- \n")
    print(" Identifier -> (accession number / GI number) \n Database -> (nucleotide,nuccore) \n Output ->  (fasta/genbank) \n Display ->  (text,xml)")

    userEmail = input("PLEASE ENTER YOUR E-MAIL ID ---> ")
    Entrez.email = userEmail

    # USER INPUT INFORMATION FOR RETRIEVAL SYSTEm

    print("NOW YOU ARE GOING TO ENTER RESPECTIVE INFORMATION FOR DATA ACCESS / RETRIEVAL : \n")
    # id
    Id = input("ID OF YOUR DATA [AS PER NCBI] : ")
    # db
    dbType = input(
        "WHAT IS THE TYPE OF DATABASE YOU WANT TO ACCESS [nucleotide / nuccore] : ")
    # rettype
    outType = input("OUTPUT FILE FORMAT [fasta / genbank] : ")
    # retmode
    dispType = input("DISPLAY FILE TYPE [text / xml]")

    optSel = input(
        "DO YOU WANT TO JUST READ THE FILE OR DOWNLOAD IT -->  [R : READ || D : DOWNLOAD]")

    if optSel == "R":
        actionHandle = Entrez.efetch(
            db=dbType, id=Id, rettype=outType, retmode=dispType)
        print(actionHandle.read())
    elif optSel == "D":
        fileAddress = input(
            "ENTER THE RELATIVE PATH / DIRECTORY WHERE YOU WANT TO STORE YOUR FILE : (EXAMPLE : C:/Users/Files/filename.format)\n")
        actionHandle = Entrez.efetch(
            db=dbType, id=Id, rettype=outType, retmode=dispType)
        fileRecord = SeqIO.read(actionHandle, outType)

        fileAddress_Final = fileAddress+"/file.gb"

        # DOWNLOAD THE FILE
        SeqIO.write(fileRecord, fileAddress_Final, outType)


# KEY MAPPING
while(True):
    selectMode = input("ENTER YOUR PREFERRED ACTION ---> ").upper()
    if selectMode == "A":
        DNA_Analysis()

    elif selectMode == "B":
        Codon_Table()

    elif selectMode == "C":
        Custom_Data()

    elif selectMode == "D":
        Access_Files()

    elif selectMode == "E":
        GC_Content()

    elif selectMode == "F":
        Amino_Acid_Map()

    elif selectMode == "G":
        NCBI_Entrez()

    elif selectMode == "Q":
        break
