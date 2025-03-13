import matplotlib.pyplot as plt
import numpy as np
import json
import seaborn as sns
import os

def gen_graphics(mode, n, tol, output):

    folder = "results/" + mode + "/" + "smwtp" + "/" + str(n) + "/" + tol
    folder_arp = "results/" + mode + "/" + "arp" + "/" + str(n) + "/" + tol

    datos = []

    for file in os.listdir(folder):
        with open(folder + "/" + file) as f:
            datos.append(json.load(f))

    orders = np.linspace(0, n, n + 1)[0:-1]

    # Hacer un gráfico de boxplot para los NMAES, NPGOS y Normalized GOs

    NMAES = []
    NPGOs = []
    NGOs = []

    for order in orders:
        nmaes = []
        npgos = []
        ngos = []
        for d in datos:
            order = int(order)
            if str(order) in d["orders"]:
                nmaes.append(d["orders"][str(order)]["Normalized MAE"])
                npgos.append(d["orders"][str(order)]["NPGO"])
                ngos.append(d["orders"][str(order)]["Normalized distance GOs"])
        NMAES.append(nmaes)
        NPGOs.append(npgos)
        NGOs.append(ngos)
    
    datos = []

    for file in os.listdir(folder_arp):
        with open(folder_arp + "/" + file) as f:
            datos.append(json.load(f))

    orders = np.linspace(0, n, n + 1)[0:-1]

    # Hacer un gráfico de boxplot para los NMAES, NPGOS y Normalized GOs

    NMAES_arp = []
    NPGOs_arp = []
    NGOs_arp = []

    for order in orders:
        nmaes = []
        npgos = []
        ngos = []
        for d in datos:
            order = int(order)
            if str(order) in d["orders"]:
                nmaes.append(d["orders"][str(order)]["Normalized MAE"])
                npgos.append(d["orders"][str(order)]["NPGO"])
                ngos.append(d["orders"][str(order)]["Normalized distance GOs"])
        NMAES_arp.append(nmaes)
        NPGOs_arp.append(npgos)
        NGOs_arp.append(ngos)

    # Graficar boxplots de NMAE en escala logarítmica

    # Configuración del estilo de las gráficas
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_context("notebook", font_scale=1.2)

    # Graficar boxplots de NMAE en escala logarítmica con colores y líneas
    plt.figure(figsize=(10, 6))
    box = plt.boxplot(NMAES, positions=orders, patch_artist=True)
    box_arp = plt.boxplot(NMAES_arp, positions=orders, patch_artist=True)

    # Colores para los boxplots
    color_smwtp = 'skyblue'
    color_arp = 'lightgreen'
    for patch in box['boxes']:
        patch.set_facecolor(color_smwtp)
    for patch in box_arp['boxes']:
        patch.set_facecolor(color_arp)

    # Añadir líneas que unan los puntos
    mean_nmaes = [np.mean(nmae) for nmae in NMAES]
    plt.plot(orders, mean_nmaes, 'o-', color='blue', linewidth=2, markersize=8, label='SMWTP')

    mean_nmaes_arp = [np.mean(nmae) for nmae in NMAES_arp]
    plt.plot(orders, mean_nmaes_arp, 'o-', color='green', linewidth=2, markersize=8, label='ARP')

    plt.xlabel('Orden', fontsize=14)
    plt.ylabel('Normalized MAE' + (" (log scale)" if tol != "r" else ""), fontsize=14)
    plt.title('Boxplot de Normalized MAE vs Orden', fontsize=16)
    if tol != "r":
        plt.yscale('log')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.xlim(-0.5, n + 0.5)  # Ajustar los límites del eje x
    plt.tight_layout()
    plt.legend()
    plt.savefig(output + '/nmae_vs_order.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Graficar boxplots de NPGO con colores y líneas
    plt.figure(figsize=(10, 6))
    box = plt.boxplot(NPGOs, positions=orders, patch_artist=True)
    box_arp = plt.boxplot(NPGOs_arp, positions=orders, patch_artist=True)

    # Colores para los boxplots
    for patch in box['boxes']:
        patch.set_facecolor(color_smwtp)
    for patch in box_arp['boxes']:
        patch.set_facecolor(color_arp)

    # Añadir líneas que unan los puntos
    mean_npgos = [np.mean(npgo) for npgo in NPGOs]
    plt.plot(orders, mean_npgos, 'o-', color='blue', linewidth=2, markersize=8, label='SMWTP')

    mean_npgos_arp = [np.mean(npgo) for npgo in NPGOs_arp]
    plt.plot(orders, mean_npgos_arp, 'o-', color='green', linewidth=2, markersize=8, label='ARP')

    plt.xlabel('Orden', fontsize=14)
    plt.ylabel('NPGO', fontsize=14)
    plt.title('Boxplot de NPGO vs Orden', fontsize=16)
    plt.grid(True, alpha=0.2)
    plt.xlim(-0.5, n + 0.5)  # Ajustar los límites del eje x
    plt.ylim(0, 1.1)
    plt.tight_layout()
    plt.legend()
    plt.savefig(output + '/npgo_vs_order.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Graficar boxplots de Normalized distance GOs con colores y líneas
    plt.figure(figsize=(10, 6))
    box = plt.boxplot(NGOs, positions=orders, patch_artist=True)
    box_arp = plt.boxplot(NGOs_arp, positions=orders, patch_artist=True)

    # Colores para los boxplots
    for patch in box['boxes']:
        patch.set_facecolor(color_smwtp)
    for patch in box_arp['boxes']:
        patch.set_facecolor(color_arp)

    # Añadir líneas que unan los puntos
    mean_ngos = [np.mean(ngo) for ngo in NGOs]
    plt.plot(orders, mean_ngos, 'o-', color='blue', linewidth=2, markersize=8, label='SMWTP')

    mean_ngos_arp = [np.mean(ngo) for ngo in NGOs_arp]
    plt.plot(orders, mean_ngos_arp, 'o-', color='green', linewidth=2, markersize=8, label='ARP')

    plt.xlabel('Orden', fontsize=14)
    plt.ylabel('Normalized distance GOs', fontsize=14)
    plt.title('Boxplot de Normalized distance GOs vs Orden', fontsize=16)
    plt.grid(True, alpha=0.2)
    plt.xlim(-0.5, n + 0.5)  # Ajustar los límites del eje x
    plt.ylim(0, 1.1)
    plt.tight_layout()
    plt.legend()
    plt.savefig(output + '/norm_distance_vs_order.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    from argparse import ArgumentParser, RawDescriptionHelpFormatter, _StoreTrueAction, ArgumentDefaultsHelpFormatter, Action
    parser = ArgumentParser(description="Permutation surrogates")
    parser.add_argument('--n', type=int, help='instance file')
    parser.add_argument('--tol', type=str, help='tolerance')
    parser.add_argument("--mode", type=str, default="YSR", help="mode")
    args = parser.parse_args()
    n = args.n
    mode = args.mode
    tol = args.tol
    gen_graphics(mode, n, tol, "graphs/" + mode + "/" + str(n) + "/" + tol)


