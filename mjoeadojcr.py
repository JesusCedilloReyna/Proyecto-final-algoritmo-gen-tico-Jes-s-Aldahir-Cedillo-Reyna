import random
import pandas as pd
import blosum
import copy
import time

import matplotlib.pyplot as plt


blosum62 = blosum.BLOSUM(62)
NFE = 0
start_time = time.time()

def get_sequences():
    seq1 = "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDLYDDDDKDRWGKLVVLGAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQV"
    seq2 = "MKTLLVAAAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQKELQKQLGQKAKEL"
    seq3 = "MAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKALCVFAIN"
    return [list(seq1), list(seq2), list(seq3)]

def crear_individuo():
    return get_sequences()

def crear_poblacion_inicial(n=10):
    individuo_base = crear_individuo()
    poblacion = [ [row[:] for row in individuo_base] for _ in range(n) ]
    return poblacion

def mutar_poblacion_v2(poblacion, num_gaps=1):
    poblacion_mutada = []
    for individuo in poblacion:
        nuevo_individuo = []
        for fila in individuo:
            fila_mutada = fila[:]
            posiciones = set()
            for _ in range(num_gaps):
                pos = random.randint(0, len(fila_mutada))
                while pos in posiciones:
                    pos = random.randint(0, len(fila_mutada))
                posiciones.add(pos)
                fila_mutada.insert(pos, '-')
            nuevo_individuo.append(fila_mutada)
        poblacion_mutada.append(nuevo_individuo)
    return poblacion_mutada

def igualar_longitud_secuencias(individuo, gap='-'):
    max_len = max(len(fila) for fila in individuo)
    individuo_igualado = [fila + [gap]*(max_len - len(fila)) for fila in individuo]
    return individuo_igualado


def evaluar_individuo_blosum62(individuo):
    global NFE
    NFE += 1
    score = 0
    n_seqs = len(individuo)
    seq_len = len(individuo[0])
    for col in range(seq_len):
        for i in range(n_seqs):
            for j in range(i+1, n_seqs):
                a = individuo[i][col]
                b = individuo[j][col]
                if a == '-' or b == '-':
                    score -= 4
                else:
                    score += blosum62[a][b]
                    
    return score

def eliminar_peores(poblacion, scores, porcentaje=0.5):
    idx_ordenados = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    n_seleccionados = int(len(poblacion) * porcentaje)
    
    ind_seleccionados = [poblacion[i] for i in idx_ordenados[:n_seleccionados]]
    scores_seleccionados = [scores[i] for i in idx_ordenados[:n_seleccionados]]
    
    return ind_seleccionados, scores_seleccionados

def cruzar_individuos_doble_punto(ind1, ind2):
    hijo1 = []
    hijo2 = []
    for seq1, seq2 in zip(ind1, ind2):
        aa_indices = [i for i, a in enumerate(seq1) if a != '-']
        if len(aa_indices) < 6:
            hijo1.append(seq1[:])
            hijo2.append(seq2[:])
            continue

        # Asegurar que el segmento tenga al menos 5 aminoácidos
        intentos = 0
        while True:
            p1, p2 = sorted(random.sample(aa_indices, 2))
            if p2 - p1 >= 5 or intentos > 10:
                break
            intentos += 1

        def cruza(seqA, seqB):
            aaA = [a for a in seqA if a != '-']
            aaB = [a for a in seqB if a != '-']
            nueva = aaA[:p1] + aaB[p1:p2] + aaA[p2:]
            resultado = []
            idx = 0
            for a in seqA:
                if a == '-':
                    resultado.append('-')
                else:
                    resultado.append(nueva[idx])
                    idx += 1
            return resultado

        nueva_seq1 = cruza(seq1, seq2)
        nueva_seq2 = cruza(seq2, seq1)

        hijo1.append(nueva_seq1)
        hijo2.append(nueva_seq2)

    hijo1 = mutar_individuo(hijo1, 1, 0.8)
    hijo2 = mutar_individuo(hijo2, 1, 0.8)
    return hijo1, hijo2


def mutar_individuo(individuo, n_gaps, p):
    """
    Mutar un solo individuo insertando n_gaps en posiciones aleatorias de sus secuencias con probabilidad p.
    individuo: lista de listas (secuencias)
    n_gaps: número de gaps a insertar por secuencia
    p: probabilidad de insertar gap en cada secuencia
    """
    nuevo_individuo = []
    for secuencia in individuo:
        sec = secuencia[:]
        if random.random() < p:
            posiciones = set()
            for _ in range(n_gaps):
                pos = random.randint(0, len(sec))
                while pos in posiciones:
                    pos = random.randint(0, len(sec))
                posiciones.add(pos)
                sec.insert(pos, '-')
        nuevo_individuo.append(sec)
    return nuevo_individuo

def cruzar_poblacion_doble_punto(poblacion):
    """
    Selecciona parejas aleatorias y realiza cruza de doble punto de corte ignorando gaps.
    Duplica el tamaño de la población.
    """
    nueva_poblacion = []
    n = len(poblacion)
    indices = list(range(n))
    random.shuffle(indices)
    parejas = [(indices[i], indices[i+1]) for i in range(0, n-1, 2)]
    if n % 2 == 1:
        parejas.append((indices[-1], indices[0]))
    for idx1, idx2 in parejas:
        padre1 = poblacion[idx1]
        padre2 = poblacion[idx2]
        hijo1, hijo2 = cruzar_individuos_doble_punto(padre1, padre2)
        
#         score_padre1 = evaluar_individuo_blosum62(padre1)
#         score_padre2 = evaluar_individuo_blosum62(padre2)
#         score_hijo1 = evaluar_individuo_blosum62(hijo1)
#         score_hijo2 = evaluar_individuo_blosum62(hijo2)

        #print(f"Padre1: {score_padre1}, Padre2: {score_padre2}, Hijo1: {score_hijo1}, Hijo2: {score_hijo2}")
        
        nueva_poblacion.append(copy.deepcopy(padre1))
        nueva_poblacion.append(copy.deepcopy(padre2))
        nueva_poblacion.append(hijo1)
        nueva_poblacion.append(hijo2)
    # Recortar si se excede el doble
    return nueva_poblacion[:2*n]

def validar_poblacion_sin_gaps(poblacion, originales):
    """
    Valida que, al eliminar los gaps de cada secuencia de cada individuo,
    las secuencias resultantes sean idénticas a las originales.
    Devuelve True si todas coinciden, False si alguna no.
    """
    for individuo in poblacion:
        for seq, seq_orig in zip(individuo, originales):
            seq_sin_gaps = [a for a in seq if a != '-']
            seq_orig_sin_gaps = [a for a in seq_orig if a != '-']
            if seq_sin_gaps != seq_orig_sin_gaps:
                return False
    return True





def obtener_best(scores, poblacion):

    if len(scores) == 0:
        return None, float('-inf')

    idx_mejor = scores.index(max(scores))

    fitness_best = scores[idx_mejor]

    best = copy.deepcopy(poblacion[idx_mejor])

    return best, fitness_best


if __name__ == "__main__":
    veryBest = None
    fitnessVeryBest = None
    poblacion = crear_poblacion_inicial(10)
    poblacion = mutar_poblacion_v2(poblacion, num_gaps=1)
    poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]
    scores = [evaluar_individuo_blosum62(ind) for ind in poblacion]
    poblacion, scores = eliminar_peores(poblacion, scores)
    
    
    for generaciones in range(100):
        poblacion = cruzar_poblacion_doble_punto(poblacion)
        #poblacion = mutar_poblacion_v2(poblacion, num_gaps=1)
        poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]
        scores = [evaluar_individuo_blosum62(ind) for ind in poblacion]
        #print(scores)
        poblacion, scores = eliminar_peores(poblacion, scores)
        #print(scores)
        best, fitness_best = obtener_best(scores, poblacion)
        if veryBest is None or fitness_best>fitnessVeryBest:
            veryBest = best
            fitnessVeryBest = fitness_best
        end_time = time.time()
        transcurrido = end_time - start_time
        print("fitness: ", fitnessVeryBest, "NFE: ", NFE, "time: ", transcurrido )
    print(veryBest)
    print ("validacion de integridad: ",validar_poblacion_sin_gaps(poblacion, get_sequences()))
    
    
def seleccion_torneo(poblacion, scores, k=3):

    if len(poblacion) == 0:
        return []

    k = min(k, len(poblacion))

    seleccionados = []

    for _ in range(len(poblacion)):

        candidatos = random.sample(
            list(zip(poblacion, scores)),
            k
        )

        ganador = max(candidatos, key=lambda x: x[1])

        seleccionados.append(
            copy.deepcopy(ganador[0])
        )

    return seleccionados

def mutacion_hibrida(individuo, p_insert=0.3, p_mover=0.3, p_eliminar=0.3):
    nuevo = []
    for sec in individuo:
        sec = sec[:]
        # Insertar
        if random.random() < p_insert:
            pos = random.randint(0, len(sec))
            sec.insert(pos, '-')
        # Mover
        if '-' in sec and random.random() < p_mover:
            idx = sec.index('-')
            sec.pop(idx)
            nueva_pos = random.randint(0, len(sec))
            sec.insert(nueva_pos, '-')
        # Eliminar
        if '-' in sec and random.random() < p_eliminar:
            sec.remove('-')
        nuevo.append(sec)
    return nuevo

def cruzar_individuos_multipunto(ind1, ind2, cortes=3):
    hijo1, hijo2 = [], []
    for s1, s2 in zip(ind1, ind2):
        aa1 = [a for a in s1 if a != '-']
        aa2 = [a for a in s2 if a != '-']
        
        if len(aa1) <= cortes + 1:
            hijo1.append(s1[:])
            hijo2.append(s2[:])
            continue
        puntos = sorted(random.sample(range(1, len(aa1)-1), cortes))
        segmentos = []
        prev = 0
        for p in puntos + [len(aa1)]:
            segmentos.append(aa1[prev:p])
            prev = p
        # alternar segmentos
        nueva1, nueva2 = [], []
        for i, seg in enumerate(segmentos):
            if i % 2 == 0:
                nueva1 += seg
                nueva2 += aa2[sum(len(s) for s in segmentos[:i]):sum(len(s) for s in segmentos[:i+1])]
            else:
                nueva1 += aa2[sum(len(s) for s in segmentos[:i]):sum(len(s) for s in segmentos[:i+1])]
                nueva2 += seg
        # reconstruir con gaps
        def reconstruir(base, nueva):
            res, idx = [], 0
            for a in base:
                if a == '-':
                    res.append('-')
                else:
                    res.append(nueva[idx]); idx += 1
            return res
        hijo1.append(reconstruir(s1, nueva1))
        hijo2.append(reconstruir(s2, nueva2))
    return hijo1, hijo2

def algoritmo_mejorado(generaciones=100, poblacion_inicial=10):
    originales = get_sequences()
    poblacion = crear_poblacion_inicial(poblacion_inicial)
    poblacion = mutar_poblacion_v2(poblacion, num_gaps=1)
    poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]
    scores = [evaluar_individuo_blosum62(ind) for ind in poblacion]
    best_hist = []

    veryBest, fitnessVeryBest = None, float('-inf')

    for gen in range(generaciones):
        # Elitismo
        elite, elite_score = obtener_best(scores, poblacion)

        # Selección torneo
        poblacion = seleccion_torneo(poblacion, scores)

        # Cruza multipunto
        nueva_poblacion = []
        for i in range(0, len(poblacion)-1, 2):
            h1, h2 = cruzar_individuos_multipunto(poblacion[i], poblacion[i+1])
            nueva_poblacion += [h1, h2]
        poblacion += nueva_poblacion

        # Mutación híbrida
        poblacion = [mutacion_hibrida(ind) for ind in poblacion]

        # Igualar longitudes
        poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]


        # Reintroducción cada 20 generaciones
        if gen % 20 == 0:
            nuevos= crear_poblacion_inicial(2)
            nuevos = mutar_poblacion_v2(nuevos, num_gaps=1)
            nuevos = [igualar_longitud_secuencias(ind) for ind in nuevos]
            poblacion += nuevos
            
        # Evaluar
        scores = [evaluar_individuo_blosum62(ind) for ind in poblacion]
        poblacion, scores = eliminar_peores(poblacion, scores)


        # Elitismo: conservar mejor
        peor_idx= scores.index(min(scores))
        poblacion[peor_idx]= elite
        scores[peor_idx]= elite_score

        best, fitness_best = obtener_best(scores, poblacion)
        if fitness_best > fitnessVeryBest:
            veryBest, fitnessVeryBest = best, fitness_best
        best_hist.append(fitness_best)

        print(f"Gen {gen}: Mejor fitness {fitness_best}")

    print("Validación de integridad:", validar_poblacion_sin_gaps(poblacion, originales))
    return best_hist


hist_original = []  
for generaciones in range(100):
    poblacion = cruzar_poblacion_doble_punto(poblacion)
    poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]
    scores = [evaluar_individuo_blosum62(ind) for ind in poblacion]
    poblacion, scores = eliminar_peores(poblacion, scores)
    best, fitness_best = obtener_best(scores, poblacion)
    if veryBest is None or fitness_best > fitnessVeryBest:
        veryBest = best
        fitnessVeryBest = fitness_best
    hist_original.append(fitness_best)
hist_mejorado = algoritmo_mejorado()

plt.plot(hist_original, label="Original")
plt.plot(hist_mejorado, label="Mejorado")
plt.xlabel("Generación")
plt.ylabel("Mejor Fitness")
plt.title("Comparación Algoritmo Genético Original vs Mejorado")
plt.legend()
plt.show()

print("\nMejor fitness ORIGINAL:", max(hist_original)) 
print("Mejor fitness MEJORADO:", max(hist_mejorado))


    
  
    
   


