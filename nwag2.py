from operator import itemgetter

class NeedlemanWunschAffineGap:

  def __init__():
    UP = (-1,0)
    LEFT = (0, -1)
    TOPLEFT = (-1, -1)
    GAPHOP = (0, 0)
    ORIGIN = (0, 0)
    D_TAG = "D"
    LG_TAG = "LG"
    UG_TAG = "UG"

  def traceback_global(v, w, pointers):
      i,j = len(v), len(w)
      curr_tracker = D_TAG
      new_v = []
      new_w = []
      while True:
          di, dj , curr_tracker = pointers[curr_tracker][i][j]
          if (di,dj) == GAPHOP:
              continue
          elif (di,dj) == LEFT:
              new_v.append('-')
              new_w.append(w[j-1])
          elif (di,dj) == UP:
              new_v.append(v[i-1])
              new_w.append('-')
              curr_tracker = UG_TAG
          elif (di,dj) == TOPLEFT:
              new_v.append(v[i-1])
              new_w.append(w[j-1])
              curr_tracker = D_TAG
  
          i, j = i + di, j + dj
          if (i <= 0 and j <= 0):
              break
      return ''.join(new_v[::-1])+'\n'+''.join(new_w[::-1])
  
  
  #phi is the gap start penalty
  #sigma is the gap continuation penalty
  def global_align(v, w, delta, phi, sigma):
      """
      Returns the score of the maximum scoring alignment of the strings v and w, as well as the actual alignment as
      computed by traceback_global.
  
      :param: v
      :param: w
      :param: delta
      """
      #Left gap
      LG = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
      #Diagnal
      D = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
      #Up gap
      UG = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
  
  
      gap_trackers = [LG_TAG, D_TAG, UG_TAG]
      pointers = {gt:[[ORIGIN for j in range(len(w)+1)] for i in range(len(v)+1)] for gt in gap_trackers}
      score, alignment = None, None
  
      d2u = GAPHOP + (UG_TAG,)
      d2l = GAPHOP + (LG_TAG,)
      d2d = TOPLEFT + (D_TAG,)
      l2l = LEFT + (LG_TAG,)
      l2d = LEFT + (D_TAG,)
      u2u = UP + (UG_TAG,)
      u2d = UP + (D_TAG,)
  
      pointer_list = {}
      pointer_list[D_TAG] = [d2u, d2l, d2d]
      pointer_list[LG_TAG] = [l2l, l2d]
      pointer_list[UG_TAG] = [u2u, u2d]
  
      #Initialize LG
      for i in range(1, len(LG[0])):
        LG[0][i] = -1* (phi + (sigma*i))
        pointers[LG_TAG][0][i] = l2l
  
        UG[0][i] = float('-inf')
      #Initialize UG
      for i in range(1, len(UG)):
        UG[i][0] = -1* (phi + (sigma*i))
        pointers[UG_TAG][i][0] = u2u
  
        LG[i][0] = float('-inf')
      #Initialize D
      D[0][0] == 0
      for i in range(1, len(D[0])):
        D[0][i] = LG[0][i]
        pointers[D_TAG][0][i] = d2l
      for i in range(1, len(D)):
        D[i][0] = UG[i][0]
        pointers[D_TAG][i][0] = d2u
  
      def get_max_step(pointers, scores):
        max_pointer, max_score = max( zip(pointers, scores), key=itemgetter(1))
        return max_pointer, max_score
  
      for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
          #Compute max for LG
          l = LG[i][j-1] - sigma
          d = D[i][j-1] - (phi + sigma)
          scores = [l, d]
          pointers[LG_TAG][i][j], LG[i][j] = get_max_step(pointer_list[LG_TAG], scores)
          #Compute max for UG
          u = UG[i-1][j] - sigma
          d = D[i-1][j] - (phi + sigma)
          scores = [u, d]
          pointers[UG_TAG][i][j], UG[i][j] = get_max_step(pointer_list[UG_TAG], scores)
          #Compute max for D
          u = UG[i][j]
          l = LG[i][j]
          d = D[i-1][j-1] + delta[v[i-1]][w[j-1]]
          scores = [u, l, d]
          pointers[D_TAG][i][j], D[i][j] = get_max_step(pointer_list[D_TAG], scores)
      score = D[-1][-1]
  
      alignment = traceback_global(v,w, pointers)
      return score, alignment
