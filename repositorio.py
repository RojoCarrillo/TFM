import csv

# A9pertura de lncRNA con todos los transcritos:
with open('lncRNA_con_transcritos.csv', newline='') as genes_lncRNA_posiciones:
    genes_lncRNA_posiciones_reader = csv.reader(genes_lncRNA_posiciones)
    nombres_lncRNA_biomart=list()
    distancias={}
    posicion_as={}
    genes_y_transcritos={}
    antisentidos=[]
    cadena_forward=[]
    diccionario_para_tabla={}
    cadena_reverse=[]
    posicion_as_cromosoma={}

    for row in genes_lncRNA_posiciones_reader: # Cada fila contiene los datos que se recogen en variables a continuación.
        nombre = row[6] #Columna con los nombres
        posicion_inicio = row[2]
        posicion_final = row[3]
        cadena = row[4]
        cromosoma= row[1]
        id=row[0]
        transcrito = row[8]
        if len(nombre) > 0: #Seleccionar los genes que tienen como producto un ARNlnc 
            nombres_lncRNA_biomart.append(nombre)
        if len(cromosoma) < 3: #Seleccionar los genes que tienen posición cromosómica definida. 
            for i in range(len(nombre)):
                nombres=nombre[i:i+3] #Rango de la longitud del nombre para poder selccionar los genes que contienen -AS en el nombre
                if nombres == "-AS" and nombre in genes_y_transcritos: 
                    genes_y_transcritos[nombre].append(transcrito)
                elif nombres == "-AS" and nombre not in genes_y_transcritos:
                    genes_y_transcritos[nombre]=[transcrito]
                if nombres == "-AS" and nombre not in antisentidos:  #Si contiene -AS en el nombre añadirlo a la lista que contiene todos los -AS
                    tamaño = int(posicion_final) - int(posicion_inicio) #Del antisentido seleccionado crear una variable con su tamaño
                    distancias[nombre]= tamaño #Y añadirlo al diccionario que recoge los tamaños de cada gen para comapararlos posteriormente
                    posicion_as[nombre]=[posicion_inicio,posicion_final]
                    posicion_as_cromosoma[nombre]=[cromosoma,posicion_inicio,posicion_final,nombre]
                    antisentidos.append(nombre)
                    diccionario_para_tabla[nombre]=[cromosoma,posicion_inicio,posicion_final,cadena,tamaño,id,nombre] #Se crea este diccionario para recoger los datos y posteriormente poder imprimirlos en formato tsv para su manejo con bash 
                    if cadena == "1" and nombre not in cadena_forward:
                        cadena_forward.append(nombre)
                    elif cadena == "-1" and nombre not in cadena_reverse:
                        cadena_reverse.append(nombre)


for gen in genes_y_transcritos:
    transcrito=genes_y_transcritos[gen] #Se lee el diccionario donde se han añadido las veces que aparece un gen para saber su número de transcritos
    if gen in diccionario_para_tabla:
        diccionario_para_tabla[gen].append(len(transcrito)) 


for gen in diccionario_para_tabla: #Los datos filtrados anteriormente recogidos en este diccionario se imprimen en formato tsv para su posterior manejo.
    chr=diccionario_para_tabla[gen][0]
    pos_init=diccionario_para_tabla[gen][1]
    pos_fin=diccionario_para_tabla[gen][2]
    strand=diccionario_para_tabla[gen][3]
    size=diccionario_para_tabla[gen][4]
    id=diccionario_para_tabla[gen][5]
    name=diccionario_para_tabla[gen][6]
    number_trans=diccionario_para_tabla[gen][7]

    #print(chr,pos_init,pos_fin,strand,size,id,name,number_trans,sep="\t")

antisentido=list() # Inicializar la variable donde se vuelven a anotar los genes antisentido con un filtrado.

# Genes AS según su nomenclatura numérica
antisentido_gene_as_total=list()
antisentido_gene_as1_total=list()
antisentido_gene_as2_total=list()
antisentido_gene_as3_total=list()
antisentido_gene_as4_total=list()
antisentido_gene_as5_total=list()
antisentido_gene_as6_total=list()
antisentido_gene_as7_total=list()


for nombre in antisentidos:
    for i in range(len(nombre)):
        nombres=nombre[i:i+3]
        if nombres == "-AS" and nombre not in antisentido:
            antisentido.append(nombre)
            

for nombre in antisentido:
    for i in range(len(nombre)):
        nombres=nombre[i:i+4]    
        if nombres == "-AS7":
            antisentido_gene_as7_total.append(nombre)

for nombre in antisentido:
    for i in range(len(nombre)):
        nombres=nombre[i:i+4]   
        if nombres == "-AS6":
            antisentido_gene_as6_total.append(nombre)                

for nombre in antisentido:
    for i in range(len(nombre)):
        nombres=nombre[i:i+4]        
        if nombres == "-AS5":
            antisentido_gene_as5_total.append(nombre)                

for nombre in antisentido:
    for i in range(len(nombre)):
        nombres=nombre[i:i+4]                     
        if nombres == "-AS4":
            antisentido_gene_as4_total.append(nombre)                
    
for nombre in antisentido:
    for i in range(len(nombre)):
        nombres=nombre[i:i+4]                     
        if nombres == "-AS3":
            antisentido_gene_as3_total.append(nombre)                
                        
for nombre in antisentido:
    for i in range(len(nombre)):
        nombres=nombre[i:i+4]
        if nombres == "-AS2":
            antisentido_gene_as2_total.append(nombre)

for nombre in antisentido:
    for i in range(len(nombre)):
        nombres=nombre[i:i+4]
        if nombres == "-AS1":
            antisentido_gene_as1_total.append(nombre)


for nombre in antisentido:
    for i in range(len(nombre)):
        nombres=nombre[i:i+3]                     
        if nombres== "-AS" and nombre not in antisentido_gene_as1_total and nombre not in antisentido_gene_as2_total and nombre not in antisentido_gene_as3_total and nombre not in antisentido_gene_as4_total and nombre not in antisentido_gene_as5_total and nombre not in antisentido_gene_as6_total:
            antisentido_gene_as_total.append(nombre)


print("La cantidad de antisentidos es", len(antisentido))
print("Los genes que tienen antisentido -AS son",len(antisentido_gene_as_total))
print("Genes que contienen -AS1 en total",len(antisentido_gene_as1_total))
print("Genes que contienen -AS2 en total",len(antisentido_gene_as2_total))
print("Genes que contienen -AS3 en total",len(antisentido_gene_as3_total))
print("Genes que contienen -AS4 en total",len(antisentido_gene_as4_total))
print("Los genes -AS4 son",sorted(antisentido_gene_as4_total))
print("Genes que contienen -AS5 en total",len(antisentido_gene_as5_total))
print("Los genes -AS5 son",sorted(antisentido_gene_as5_total))
print("Genes que contienen -AS6 en total",len(antisentido_gene_as6_total))
print("Los genes -AS6 son",sorted(antisentido_gene_as6_total))
print("La cantidad de antisentidos -AS7 es", len(antisentido_gene_as7_total))
print("----------------------------------------------------------------")



genes_con_as_nueva=[]

#Listas con los genes sentido teoricos
genes_AS=[]
genes_AS1=[]
genes_AS2=[]
genes_AS3=[]
genes_AS4=[]
genes_AS5=[]
genes_AS6=[]

comprobacion=[]

as_1=[] #Inicialización de dos listas donde se recogerá su gen sentido teorico según la definición
as_2=[]

for gen in antisentido_gene_as_total:
    gen1=gen[0:(len(gen)-3)] #Se elimina del gen AS la parte correspondiente a -AS
    genes_con_as_nueva.append(gen1)
    genes_AS.append(gen1)
    as_1.append(gen)

for gen in antisentido_gene_as1_total:
    gen1=gen[0:(len(gen)-4)]
    genes_con_as_nueva.append(gen1)
    genes_AS1.append(gen1)
    as_2.append(gen)

for gen in antisentido_gene_as2_total:
    gen1=gen[0:(len(gen)-4)]
    genes_AS2.append(gen1)
    as_2.append(gen)

for gen in antisentido_gene_as3_total:
    gen1=gen[0:(len(gen)-4)]
    genes_con_as_nueva.append(gen1)
    genes_AS3.append(gen1)
    as_2.append(gen)

for gen in antisentido_gene_as4_total:
    gen1=gen[0:(len(gen)-4)]
    genes_con_as_nueva.append(gen1)
    genes_AS4.append(gen1)
    as_2.append(gen)

for gen in antisentido_gene_as5_total:
    gen1=gen[0:(len(gen)-4)]
    genes_con_as_nueva.append(gen1)
    genes_AS5.append(gen1)
    as_2.append(gen)

for gen in antisentido_gene_as6_total:
    gen1=gen[0:(len(gen)-4)]
    genes_con_as_nueva.append(gen1)
    genes_AS6.append(gen1)
    as_2.append(gen)

filtro=[]
for gen in genes_con_as_nueva:
    if gen not in filtro:
        filtro.append(gen)

with open('protein_coding_con_transcritos.csv', newline='') as genes_protein_coding_posiciones: #Mismo procedimiento que el caso anterior de lectura de ficheros para los genes codificantes de proteínas
    genes_protein_coding_posiciones_reader = csv.reader(genes_protein_coding_posiciones)
    nombres_protein_coding_biomart=list()
    genes_con_as=list()
    distancias_genes={}
    posicion_genes={}
    cadena_forward_gen=[]
    cadena_reverse_gen=[]
    genes_y_transcritos_prot_coding={}
    genes_y_transcritos_prot_coding_con_as={}
    posicion_genes_con_as={}
    control=[]
    tabla_final_genes_prot_coding={}
    for row in genes_protein_coding_posiciones_reader:
        nombre = row[6]
        posicion_inicio = row[2]
        posicion_final = row[3]
        cadena = row[4]
        id=row[0]
        cromosoma= row[1]
        transcrito = row[8]
        if len(cromosoma) < 3:
            nombres_protein_coding_biomart.append(nombre)                   
            
            if nombre in genes_y_transcritos_prot_coding_con_as :
                    genes_y_transcritos_prot_coding_con_as[nombre].append(transcrito)
            elif nombre not in genes_y_transcritos_prot_coding_con_as:
                    genes_y_transcritos_prot_coding_con_as[nombre]=[transcrito] 
        
            if cadena == "1" and nombre not in cadena_forward_gen:
                cadena_forward_gen.append(nombre)
            elif cadena == "-1" and nombre not in cadena_reverse_gen:
                cadena_reverse_gen.append(nombre)

            if nombre in filtro and nombre not in tabla_final_genes_prot_coding: #Se añade a la lista los genes codificantes que tienen una coincidencia exacta con la definición de su gen antisentido. 
                distancia=int(posicion_final) - int(posicion_inicio)
                distancias_genes[nombre]=distancia
                posicion_genes[nombre]=[posicion_inicio,posicion_final]
                posicion_genes_con_as[nombre]=[cromosoma,posicion_inicio,posicion_final,nombre]
                genes_con_as.append(nombre)
                control.append(nombre)
                tabla_final_genes_prot_coding[nombre]=[cromosoma,posicion_inicio,posicion_final,cadena,distancia,id,nombre]

print(len(tabla_final_genes_prot_coding))


no_coincidencias=[] #Comprobar el desfase de numeros y porque ocurre
for gen in filtro:
    if gen not in tabla_final_genes_prot_coding:
        no_coincidencias.append(gen)

print(no_coincidencias) 
print(len(no_coincidencias))
#Aparecen nombres de genes que no codifican proteinas y otros que no se relaciona exactamente con el nombre del antisentido
#Se necesita comprobación manual para preparar el fichero de datos. 

## Tamaños

mas_grande_as=[]
mas_grande_gen=[]
iguales=[]

for clave in distancias: #Se comparan los dos diccionarios creados en la lectura de los ficheros para comprobar el tamaño de los genes sentido-antisentido con coincidencia perfecta
    for nombre in distancias_genes:
        clave1=clave[0:len(nombre)]
        if nombre == clave1:
            if distancias_genes[nombre] > distancias[clave]:
                mas_grande_gen.append(clave)
            elif distancias_genes[nombre] < distancias[clave]: 
                mas_grande_as.append(clave)
            elif distancias_genes[nombre] == distancias[clave]:
                iguales.append(clave)

print(len(mas_grande_as))
print(len(mas_grande_gen))
print(len(iguales))



