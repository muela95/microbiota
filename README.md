# microbiota
Sólo el código. Empieza siguiendo la guía de DADA2;
* Importa raws, F y R
* Filtra y trimea (está comentado lo que yo hice de primeras, sin comentar la sugerencia de Thomaz de quitar un número fijo, no solo primers)
* Aprende error rates
* Quita réplicas
* Aplica algoritmo dada
* Une pares reverse y forward
* Hace tabla de secuencias
* Quita quimeras
* Asigna taxonomía (necesitas archivo silva taxonomía)
* Alinea secuencias
* Crea árbol filogenético
* Optimiza el fitness del modelo
* Crea archivo phyloseq
* A partir de aquí ya es codigo desde la guía de Tjazi.
* Transformo ps a count tables en niveles de género y familia
* Filtro todo lo que tenga una prevalencia de menos del 10%
* Transformación CLR
* Stacked barplots
* Alfa diversidad, plots, estadística
* Beta diversidad, distancia euclídea (PCAs, por recomendación de Thomaz), permanovas
* Expresión diferencial a base de anovas individuales en cada género (suena a chapuza, intento hacer después un ALDEX2)
* Correlaciones de las bacterias en las que salen anovas significativos con conducta
