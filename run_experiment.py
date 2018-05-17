import os
import logging
import subprocess
import numpy as np
from tqdm import tqdm
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
logger = logging.getLogger(__name__)


def run_experiment(args):
    seed, logdir = args
    logpath = os.path.join(logdir, '%04d.ltsv' % (seed))
    args = [
        'java',
        'RoadsAndJunctionsVis',
        '-novis',
        '-exec',
        './a.out',
        '-seed',
        str(seed),
    ]
    with open(logpath, 'w') as fp:
        res = subprocess.run(args, stdout=fp)
    if res.returncode != 0:
        logger.error('error on {}'.format(res.args))


def parse_log(args):
    seed, logdir = args
    logpath = os.path.join(logdir, '%04d.ltsv' % (seed))

    score = 1e10
    base_score = 1e10
    with open(logpath) as fp:
        for row in fp:
            d = {}
            if row.find('Score =') >= 0:
                score = float(row.split('=')[1].strip())
                if score < 0:
                    score = 1e10
                continue
            if row.find(':') < 0:
                continue
            for chunk in row.strip().split('\t'):
                k, v = chunk.split(':')
                d[k] = v
            if 'base_score' in d:
                base_score = float(d['base_score'])

    assert score < 1e10, 'seed: {}'.format(seed)
    assert base_score < 1e10, 'seed: {}'.format(seed)

    return 100 * (base_score-score) / base_score


def run_experiments():
    now = datetime.now()
    logdir = os.path.join('./log', now.strftime('%Y%m%d-%H%M%S'))
    os.makedirs(logdir, exist_ok=True)
    executor = ProcessPoolExecutor()
    args_list = []
    for i in range(100):
        args_list.append([i+1, logdir])
    total = len(args_list)
    for res in tqdm(
        executor.map(run_experiment, args_list), total=total, ncols=0
    ):
        pass

    scores = []
    for res in tqdm(executor.map(parse_log, args_list), total=total, ncols=0):
        scores.append(res)

    scores = np.array(scores)
    mean = np.mean(scores)
    std = np.std(scores)
    max_score = np.max(scores)
    max_score_idx = np.argmax(scores) + 1
    logger.info('mean: {}, std: {}, max: {}(seed:{})'.format(mean, std, max_score, max_score_idx))


def main():
    logging.basicConfig(level=logging.INFO)
    run_experiments()


if __name__ == '__main__':
    main()
