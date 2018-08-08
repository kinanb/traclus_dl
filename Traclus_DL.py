from Trajectory import *
import sys
import datetime
from heapq import *
#import matplotlib.pyplot as plt
from collections import defaultdict
from numpy import arange
#from sets import Set
from itertools import count
step_size = 1.0
def round_to(n, precision):
    correction = 0.5 if n >= 0 else -0.5
    return int(n/precision+correction)*precision

segment_to_line_dist =  defaultdict(lambda : defaultdict(float))
segment_to_line_closest_seg =  defaultdict(lambda : defaultdict(Trajectory.TrajectorySegment))





def segments_distance(x11, y11, x12, y12, x21, y21, x22, y22, name1="", name2=""):
  """ distance between two segments in the plane:
      one segment is (x11, y11) to (x12, y12)
      the other is   (x21, y21) to (x22, y22)
  """
  if segments_intersect(x11, y11, x12, y12, x21, y21, x22, y22, name1, name2): return 1e-8
  # try each of the 4 vertices w/the other segment
  distances = []
  distances.append(point_segment_distance(x11, y11, x21, y21, x22, y22))
  distances.append(point_segment_distance(x12, y12, x21, y21, x22, y22))
  distances.append(point_segment_distance(x21, y21, x11, y11, x12, y12))
  distances.append(point_segment_distance(x22, y22, x11, y11, x12, y12))

  return min(distances)

def segments_intersect(x11, y11, x12, y12, x21, y21, x22, y22, name1, name2):
  """ whether two segments in the plane intersect:
      one segment is (x11, y11) to (x12, y12)
      the other is   (x21, y21) to (x22, y22)
  """
  dx1 = x12 - x11
  dy1 = y12 - y11
  dx2 = x22 - x21
  dy2 = y22 - y21
  delta = dx2 * dy1 - dy2 * dx1
  
  if delta == 0: return False  # parallel segments
  s = (dx1 * (y21 - y11) + dy1 * (x11 - x21)) / delta
  t = (dx2 * (y11 - y21) + dy2 * (x21 - x11)) / (-delta)
  
  return (0 <= s <= 1) and (0 <= t <= 1)

import math


def point_segment_distance(px, py, x1, y1, x2, y2):
  dx = x2 - x1
  dy = y2 - y1
  if dx == dy == 0:  # the segment's just a point
    return math.hypot(px - x1, py - y1)

  # Calculate the t that minimizes the distance.
  t = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy)

  # See if this represents one of the segment's
  # end points or a point in the middle.
  if t < 0:
    dx = px - x1
    dy = py - y1
  elif t > 1:
    dx = px - x2
    dy = py - y2
  else:
    near_x = x1 + t * dx
    near_y = y1 + t * dy
    dx = px - near_x
    dy = py - near_y

  return math.hypot(dx, dy)






def closest_point(px, py, x1, y1, x2, y2):
    
    dx = x2 - x1
    dy = y2 - y1
    if dx == dy == 0:  # the segment's just a point
        return (x1,y1)
  
  # Calculate the t that minimizes the distance.
    t = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy)

  # See if this represents one of the segment's
  # end points or a point in the middle.
    #print "closest_point:", px, py, "to: ", x1, y1, x2, y2,t,  x1+t*dx, y1+t*dy,  math.hypot(x1+t*dx-px, y1+t*dy-py)
    if t < 0:
        return (x1,y1)
    elif t > 1:
        return (x2,y2)
    else:
        near_x = x1 + t * dx
        near_y = y1 + t * dy

        return (near_x, near_y)

def sum_pairwise(segments):
    sum_pair = 0.0
    for ind1 in range(len(segments)):
        for ind2 in range(ind1+1, len(segments)):
            sum_pair += hypot(segments[ind1].startx - segments[ind2].startx, segments[ind1].starty - segments[ind2].starty)
    return sum_pair




def DBScan(seg1, traj_angles,  max_dist, min_weight, max_angle):
    #print seg1.parent.name, seg1, seg1.id, datetime.datetime.now()
#    reachable_lines = []
    reachable_segs = []
    sumweight = 0.0
    represented_lines = set()
    for angle in traj_angles:

        if abs(angle - seg1.parent.angle) > max_angle:
            continue

        for line2 in traj_angles[angle]:
            
            if(line2.name not in segment_to_line_dist[seg1.id]):
                segment_to_line_dist[seg1.id][line2.name] =  segments_distance(seg1.startx, seg1.starty, seg1.endx, seg1.endy, line2.startx, line2.starty, line2.endx, line2.endy)
                closest_x, closest_y = closest_point((seg1.startx+seg1.endx)/2.0, (seg1.starty+seg1.endy)/2.0, line2.startx, line2.starty, line2.endx, line2.endy)

                segment_to_line_closest_seg[seg1.id][line2.name] = line2.get_segment_at(closest_x, closest_y)

            if segment_to_line_dist[seg1.id][line2.name] <= max_dist:
               # print "dist_obtained: ", segments_distance(seg1.startx, seg1.starty, seg1.endx, seg1.endy, line2.startx, line2.starty, line2.endx, line2.endy) 
                represented_lines.add(line2)
  
                reachable_segs.append(segment_to_line_closest_seg[seg1.id][line2.name])
                sumweight = sumweight + line2.weight

 
    if sumweight < min_density:
        return (-1, [])
    else:
#        print "dbscan:", sumweight, seg1, reachable_segs
        return expand_cluster(seg1, traj_angles, reachable_segs, represented_lines, max_dist, min_weight, max_angle)
        

def expand_cluster(seg1, traj_angles, reachable_segs, represented_lines, max_dist, min_weight, max_angle):
    corridor_assignment = [seg1]
    sum_weight = seg1.parent.weight
    while len(reachable_segs) > 0:
        new_candidates = [] 
        for seg2 in reachable_segs: 
            
            if seg2 not in corridor_assignment:

                corridor_assignment.append(seg2)
                sum_weight += seg2.parent.weight
                sumweight_seg2 = 0.
                new_reachable = [] #for 'RegionQuery' (TODO: write separate function)
                for angle in traj_angles:
                    if abs(angle-seg1.parent.angle) > max_angle:
                        continue
                    for line3 in traj_angles[angle]:


                        if(line3.name not in segment_to_line_dist[seg2.id]):
                            segment_to_line_dist[seg2.id][line3.name] =  segments_distance(seg2.startx, seg2.starty, seg2.endx, seg2.endy, line3.startx, line3.starty, line3.endx, line3.endy)
                            closest_x, closest_y = closest_point((seg2.startx+seg2.endx)/2.0, (seg2.starty+seg2.endy)/2.0, line3.startx, line3.starty, line3.endx, line3.endy)

                            segment_to_line_closest_seg[seg2.id][line3.name] = line3.get_segment_at(closest_x, closest_y)

                        if segment_to_line_dist[seg2.id][line3.name] <= max_dist:
                            sumweight_seg2 = sumweight_seg2 + line3.weight
                            if line3 not in represented_lines:
                                represented_lines.add(line3)
                                new_reachable.append(segment_to_line_closest_seg[seg2.id][line3.name])
                if sumweight_seg2 >= min_density:
                    
                    
                    for seg3 in new_reachable:
                        if seg3 not in reachable_segs:
                            new_candidates.append(seg3)
  
        reachable_segs = new_candidates
    #print "expanded:", sum_weight, corridor_assignment, "\n"
    return (sum_weight, corridor_assignment)


pq = []                         # list of entries arranged in a heap
entry_finder = {}               # mapping of tasks to entries
REMOVED = []     # placeholder for a removed task
counter = count()     # unique sequence count
entry_hash = defaultdict(lambda : defaultdict(list))
priority_index = -1
sumsq_list = []
sumsq_index = -1
index_range = []
mindex = 0
to_removed = set()
def remove_cluster(cluster_seed):
    'Mark an existing cluster as REMOVED.  Raise KeyError if not found.'
    entry = entry_finder.pop(cluster_seed)
    entry[-1] = REMOVED



def check_removed(cluster_seed, cluster):
    removed = set()
    keeper = []
    global to_removed
    sum_weight = 0.0
    for segment in cluster:
        if segment.id in to_removed:
            #print "success with ",  segment.id
            removed.add(segment)
        else:
            keeper.append(segment)
            sum_weight += segment.parent.weight
    if len(removed):
        #print keeper, "vs: ", removed
        if sum_weight < min_density:
            remove_cluster(cluster_seed)
        else:
            sumsq = sum_pairwise(keeper)
            add_cluster(cluster_seed, keeper, int(sum_weight*100), sumsq)
        return True
    else:
        return False



def add_cluster(cluster_seed, cluster, priority, sumsq):
    'Add a new task or update the priority of an existing cluster'
    if cluster_seed in entry_finder:
        remove_cluster(cluster_seed)

    entry = [priority, sumsq, cluster_seed, cluster]
    entry_hash[priority][sumsq].append(entry)
    entry_finder[cluster_seed] = entry


def pop_cluster():
    global priority_index
    global mindex
    global seed_index
    global sumsq_index
    global sumsq_list
    if not entry_hash:
        raise KeyError('pop from an empty priority queue')
    
    if priority_index == -1:
        priority_index = max(entry_hash.keys())
        mindex = min(entry_hash.keys())
        sumsq_index = 0
        sumsq_list = sorted(entry_hash[priority_index].keys())
        seed_index = 0

    while priority_index >= mindex:
        if seed_index >= len(entry_hash[priority_index][sumsq_list[sumsq_index]]):
            seed_index = 0
            sumsq_index += 1


        if sumsq_index >= len(sumsq_list):
            if priority_index == mindex:
                raise KeyError('pop from an empty priority queue')
            priority_index -= 1
            while priority_index not in entry_hash:
                priority_index -= 1
            sumsq_list =  sorted(entry_hash[priority_index].keys())
            seed_index = 0
            sumsq_index = 0
            
       # print priority_index, sumsq_list, sumsq_index, seed_index
        priority, sumsq, cluster_seed, cluster = entry_hash[priority_index][sumsq_list[sumsq_index]][seed_index]
        if not check_removed(cluster_seed, cluster):
            for segment in cluster:
            #    print "supposed to remove", segment.id
                to_removed.add(segment.id)
            return cluster
  
        seed_index = seed_index + 1
        



def print_weighted_averages(cluster_segments, corr_number, oh):
    weightsum = 0.
    x1sum = 0.
    y1sum = 0.
    x2sum = 0.
    y2sum = 0.
    for segment in cluster_segments:
        weightsum += segment.parent.weight
        x1sum += segment.startx * segment.parent.weight
        x2sum += segment.endx * segment.parent.weight
        y1sum += segment.starty * segment.parent.weight
        y2sum += segment.endy * segment.parent.weight
    


    oh.write(str(corr_number) + "\t"+  str(weightsum) + "\tLINESTRING(" + str(x1sum/weightsum) + " "+ str(y1sum/weightsum) + ", " + str(x2sum/weightsum) + " " + str(y2sum/weightsum) + ")\n")
#main program, need to add main function

#file IO: needs to be in function
# Read file instantiate trajectories, add to dictionary of angles
infile = sys.argv[1]
max_dist = float(sys.argv[2])
min_density = float(sys.argv[3])
max_angle = float(sys.argv[4])
segment_size = float(sys.argv[5])

# For each line, get candidate lines from angle dictionary, find within distance
traj_angles = defaultdict(list) #look up angles fast, each angle (in degrees) will have a list of Trajectory objects
trajectories = [] # list of trajectories, keep it? 
fh = open (infile, 'r') 
segment_list_oh = open(infile + "." + str(max_dist) + "." + str(min_density)  + "." + str(max_angle) + "." + str(segment_size) + "_seg_list.csv", "w")
corridor_list_oh = open (infile + "." + str(max_dist) + "."  + str(min_density) + "." + str(max_angle) + "." + str(segment_size)  + "_cor_list.csv", "w")
corridor_list_oh.write("cor_id\tcor_weight\tcor_coordinates\n");
segment_list_oh.write("dl_id\tdl_weight\tdl_angle\tcor_id\tseg_coordinates\n");
for line in fh:
    linelist = line.split();
    traj = Trajectory(name=linelist[0], weight = float(linelist[1]), startx=float(linelist[2]), starty=float(linelist[3]), endx=float(linelist[4]), endy=float(linelist[5]))
    traj.make_segments(segment_size)
    rounded_angle = round_to(traj.angle, 0.5);
    traj_angles[rounded_angle].append(traj) #add the trajectory to a list in the dictionary of angles
    trajectories.append(traj)
    
    


for line in trajectories:
    for segment in line.segments:
        (sumweight, segments ) = DBScan(segment, traj_angles, max_dist, min_density, max_angle)
        if sumweight >= min_density:
            add_cluster(segment, segments, int(sumweight*100), sum_pairwise(segments))
    


corridor = 0
assigned = set()
while True:
    try:
        cluster_segments = pop_cluster()
        print_weighted_averages(cluster_segments, corridor, corridor_list_oh)


        for segment in cluster_segments:
            segment_list_oh.write(segment.parent.name + "\t" + str(segment.parent.weight) + "\t" + str(segment.parent.angle) + "\t"  + str(corridor) + "\tLINESTRING(" + str(segment.startx) + " " + str(segment.starty) + ", " + str(segment.endx) + " " + str(segment.endy) + ")\n")
            assigned.add(str(segment))
        corridor += 1
        
    except KeyError:
        break


for line in trajectories:
    for segment in line.segments:
        if str(segment) not in assigned:
            segment_list_oh.write(segment.parent.name + "\t" + str(segment.parent.weight) + "\t" + str(segment.parent.angle) + "\t"  + str(-1) + "\tLINESTRING(" + str(segment.startx) + " " + str(segment.starty) + ", " + str(segment.endx) + " " + str(segment.endy) + ")\n")




    
corridor_list_oh.close()

segment_list_oh.close()
